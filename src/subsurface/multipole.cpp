/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2010 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <boost/timer.hpp>
#include "irrtree.h"
#include "dipolebased.h"

#define NO_SSE_QUERY

MTS_NAMESPACE_BEGIN

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources. A SSE2
 * implementation doesn't make things faster here.
 */
struct IsotropicMultipoleQuery {
	inline IsotropicMultipoleQuery(const std::vector<Spectrum> &zr, const std::vector<Spectrum> &zv, 
		const Spectrum &sigmaTr, Float Fdt, const Point &p, const Normal &_ns,
		const Float _d, const Float numExtraDipoles)
		: zrList(zr), zvList(zv), sigmaTr(sigmaTr), result(0.0f), one(1.0f), d(_d), 
		  Fdt(Fdt), p(p),  ns(_ns), numExtraDipoles(numExtraDipoles), count(0) {
		// the smallest distance is at the zero order dipole
		Float zrMin = zr[numExtraDipoles + 1].min();
		zrMinSq = zrMin * zrMin;
	}

	inline void operator()(const IrradianceSample &sample) {
		const Spectrum negSigmaTr = sigmaTr * (-1.0f);
		const Spectrum rSqr = Spectrum(std::max((p - sample.p).lengthSquared(), zrMinSq));
		//const Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());

		Spectrum dMoR(0.0f);
		Spectrum dMoT(0.0f);
		for (int i=-numExtraDipoles; i<=numExtraDipoles; ++i) {
			int idx = i + numExtraDipoles;

			const Spectrum &zr = zrList[idx];
			const Spectrum &zv = zvList[idx];

			/* Distance to the real source */
			Spectrum dr = (rSqr + zr*zr).sqrt();
			/* Distance to the image point source */
			Spectrum dv = (rSqr + zv*zv).sqrt();
			Spectrum C1 = (sigmaTr + one / dr);
			Spectrum C2 = (sigmaTr + one / dv);

			/* Do not include the reduced albedo - will be canceled out later */
			dMoR += Spectrum(0.25f * INV_PI) *
					(  zr * C1 * ((negSigmaTr * dr).exp()) / (dr * dr)
					 - zv * C2 * ((negSigmaTr * dv).exp()) / (dv * dv));

			dMoT += Spectrum(0.25f * INV_PI) *
					((d - zr) * C1 * ((negSigmaTr * dr).exp()) / (dr * dr)
					 - (d + zv) * C2 * ((negSigmaTr * dv).exp()) / (dv * dv));
		}

		/* combine Mo based on R and Mo based on T to a new
		 * Mo based on a combined profile P. */
		Float cosN = dot(ns, sample.n);
		Spectrum dMoP = 0.5 * (cosN + 1) * dMoR
			+ 0.5 * (1 - cosN) * dMoT;

		result += dMoP * sample.E * (sample.area * Fdt);
	 
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	const std::vector<Spectrum> &zrList, &zvList;
	Spectrum sigmaTr, result;
	const Spectrum one;
	Spectrum d;
	Float Fdt, zrMinSq;
	Point p;
	Normal ns;
	int numExtraDipoles;
	int count;
};

/**
 * Computes the combined diffuse radiant exitance caused by a number of
 * dipole sources. It is meant to be used with a look-up-table (LUT)
 * for dMoR which is indexed with the distance. ToDo: Index with squared
 * distance to save one sqrt().
 */
struct IsotropicLUTMultipoleQuery {
	inline IsotropicLUTMultipoleQuery(const ref<LUTType> &lutR, const ref<LUTType> &lutT,
			Float _res, Float _Fdt, const Point &_p, const Normal &_ns, Float _minDist)
		: dMoR_LUT(lutR), dMoT_LUT(lutT), entries(lutR->size()), invResolution(1.0f / _res),
		  result(0.0f), Fdt(_Fdt), p(_p), ns(_ns), count(0), minDist(_minDist) {
	}

	inline void operator()(const IrradianceSample &sample) {
		//Float rSqr = (p - sample.p).lengthSquared();
		Float r = (p - sample.p).length();
		// Avoid singularities (see Jensen et al. 2001)
		r = std::max(r, minDist);
		/* Look up dMo for the distance. As in the normal query,the
		 * reduced albedo is not included. It will be canceled out
		 * later. The index is rounded to the next nearest integer. */
		int index = (int) (r * invResolution + 0.5f);
		if (index < entries) {
			Spectrum dMoR = dMoR_LUT->at(index);
			Spectrum dMoT = dMoT_LUT->at(index);

			/* combine Mo based on R and Mo based on T to a new
			 * Mo based on a combined profile P. */
			Float cosN = dot(ns, sample.n);
			Spectrum dMoP = 0.5f * (cosN + 1.0f) * dMoR
						  + 0.5f * (1.0f - cosN) * dMoT;
			result += dMoP * sample.E * (sample.area * Fdt);
		}
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	/* a reference to a dMoR look-up-table */
	const ref<LUTType> &dMoR_LUT, &dMoT_LUT;
	int entries;
	Float invResolution;
	Spectrum result;
	Float Fdt;
	Point p;
	Normal ns;
	int count;
	Float minDist;
};

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme with his and Donner's multipole approach.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *	 Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 *
 * ("Light Diffusion in Multi-Layered Tranclucent Materials", in ACM
 *	 Transactions on Graphics 24)
 */
class IsotropicMultipole : public DipoleBasedSubsurface<IsotropicMultipole> {
public:
	IsotropicMultipole(const Properties &props) 
		: DipoleBasedSubsurface(props) {
		/* The number of extra dipoles n where the total number of dipoles is
		 * 2 * n + 1. That is the dipole at the surface plus n dipoles above
		 * and below it. */
		m_extraDipoles = props.getInteger("extraDipoles", 5);
		/* The thickness of the the material slab. */
		m_slabThickness = props.getFloat("slabThickness", 1.0f);
	}
	
	IsotropicMultipole(Stream *stream, InstanceManager *manager) 
	 : DipoleBasedSubsurface(stream, manager) { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		DipoleBasedSubsurface::serialize(stream, manager);
	}

	Spectrum Mo(const Scene *scene, Sampler *sampler,  const Intersection &its,
			const Vector &d, int depth) const {
		const Normal &n = its.geoFrame.n;
		// calulate diffuse reflectance and transmittance
		IsotropicMultipoleQuery query(m_zr, m_zv, m_sigmaTr, m_Fdt, its.p, n, m_slabThickness, m_extraDipoles);
		m_octree->execute(query);

		return query.getResult();
	}

	Spectrum MoLUT(const Scene *scene, Sampler *sampler,  const Intersection &its,
			const Vector &d, int depth) const {
		const Normal &n = its.geoFrame.n;
		// calulate diffuse reflectance and transmittance
		IsotropicLUTMultipoleQuery query(m_RdLookUpTable, m_TdLookUpTable, m_lutResolution, m_Fdt, its.p, n, m_minMFP);
		m_octree->execute(query);

		return query.getResult();
	}

	void configureDipoles() {
		/* Fluence vanishing point */
		Spectrum zb = 2 * m_A * m_D;
		/* Needed for calculation purposes */
		Spectrum d = Spectrum(m_slabThickness);

		/* Distance of the dipole point sources to the surface */
		for (int i=-m_extraDipoles; i<=m_extraDipoles; ++i) {
			Float dipoleIdx = 2 * i;
			Spectrum basicDist = dipoleIdx * d + dipoleIdx * 2 * zb;
			m_zr.push_back( basicDist + m_mfp );
			m_zv.push_back( basicDist - m_mfp - 2 * zb );
		}
	}

	/* Creates a hash for the current configuraiton of this isotropic
	 * dipole. All relevant settings are included. This could e.g. be
	 * used to identify a cached look-up-table. The parameter forTd
	 * defines whether the hash is for Rd or Td.
	 */
	std::string hashConfig(bool forTd) {
		std::ostringstream oss;
		// set precision to get around floating point rounding issues
		oss.precision(5);
		oss << "multipole" << m_lutResolution << "," << m_errThreshold << "," << m_sigmaTr.toString() << ","
			<< m_alphaPrime.toString() << "," << m_extraDipoles;

		for (std::vector<Spectrum>::const_iterator it=m_zr.begin(); it!=m_zr.end(); ++it)
			oss << it->toString();
		for (std::vector<Spectrum>::const_iterator it=m_zv.begin(); it!=m_zv.end(); ++it)
			oss << it->toString();

		if (forTd)
			oss << ",d=" << m_slabThickness;

		return oss.str();
	}

	void configureLUT() {
		ref<SubsurfaceMaterialManager> smm =
			SubsurfaceMaterialManager::getInstance();
		std::string lutHashR = hashConfig(false);
		std::string lutHashT = hashConfig(true);

		if (smm->hasLUT(lutHashR) && smm->hasLUT(lutHashT)) {
			/* get Rd LUT */
			LUTRecord lutR = smm->getLUT(lutHashR);
			m_RdLookUpTable = lutR.lut;
			AssertEx(lutR.resolution == m_lutResolution,
				"Cached Rd LUT does not have requested resolution!");
			/* get Td LUT */
			LUTRecord lutT = smm->getLUT(lutHashT);
			m_TdLookUpTable = lutT.lut;
			AssertEx(lutT.resolution == m_lutResolution,
				"Cached Td LUT does not have requested resolution!");
		} else {
			// call createLUT() of base
			m_RdLookUpTable = createLUT(boost::bind(
				&mitsuba::IsotropicMultipole::getRd, this, _1));
			m_TdLookUpTable = createLUT(boost::bind(
				&mitsuba::IsotropicMultipole::getTd, this, _1));

			/* Create new LUTRecord and store this LUT */
			LUTRecord lutRecR(m_RdLookUpTable, m_lutResolution);
			smm->addLUT(lutHashR, lutRecR);
			LUTRecord lutRecT(m_RdLookUpTable, m_lutResolution);
			smm->addLUT(lutHashT, lutRecT);
			/* Make sure the LUT got actually stored. */
			AssertEx(smm->hasLUT(lutHashR), "Rd LUT is not available, but it should be!");
			AssertEx(smm->hasLUT(lutHashT), "Td LUT is not available, but it should be!");
		}
	}

	Spectrum getRd(const Spectrum r) const {
		const Spectrum one(1.0f);
		const Spectrum negSigmaTr = m_sigmaTr * (-1.0f);
		const Spectrum rSqr = r * r;

		Spectrum Rd(0.0f);
		for (int i=-m_extraDipoles; i<=m_extraDipoles; ++i) {
			int idx = i + m_extraDipoles;

			Spectrum zri = m_zr[idx];
			Spectrum zvi = m_zv[idx];
			/* Distance to the real source */
			Spectrum dri = (rSqr + zri*zri).sqrt();
			/* Distance to the image point source */
			Spectrum dvi = (rSqr + zvi*zvi).sqrt();

			/* Do not include the reduced albedo - will be canceled out later */
			Rd += Spectrum(0.25f * INV_PI) *
				 (zri * (one + m_sigmaTr * dri) * ((negSigmaTr * dri).exp())
						/ (dri * dri * dri)
				- zvi * (one + m_sigmaTr * dvi) * ((negSigmaTr * dvi).exp())
						/ (dvi * dvi * dvi));
		}
		return Rd;
	}

	Spectrum getTd(const Spectrum r) const {
		const Spectrum one(1.0f);
		const Spectrum negSigmaTr = m_sigmaTr * (-1.0f);
		const Spectrum rSqr = r * r;

		Spectrum Td(0.0f);
		for (int i=-m_extraDipoles; i<=m_extraDipoles; ++i) {
			int idx = i + m_extraDipoles;

			Spectrum zri = m_zr[idx];
			Spectrum zvi = m_zv[idx];

			/* Distance to the real source */
			Spectrum dri = (rSqr + zri*zri).sqrt();
			/* Distance to the image point source */
			Spectrum dvi = (rSqr + zvi*zvi).sqrt();
			/* Do not include the reduced albedo - will be canceled out later */
			const Spectrum d = Spectrum(m_slabThickness);
			Td += (d - zri) * (one + dri * m_sigmaTr) * ((negSigmaTr * dri).exp())
						/ (dri * dri * dri)
				- (d - zvi) * (one + dvi * m_sigmaTr) *((negSigmaTr * dvi).exp())
						/ (dvi * dvi * dvi);
		}
		/* Unfordunately, Td can get negative */
		Td.clampNegative();

		return Spectrum(0.25f * INV_PI) * Td;
	}

	MTS_DECLARE_CLASS()
private:
	/* Lists for the dipole distances */
	std::vector<Spectrum> m_zr, m_zv;
	/* The number of additional dipoles. Specifies the number of
	 * dipoles that get added below *and* above the surface. E.g.
	 * 1 additional dipole means you have the standard dipole plus
	 * one above and one below.
	 */
	int m_extraDipoles;

	Float m_slabThickness;
	/* Look-up-table for Rd and Td, indexed by the distance r. */
	ref<LUTType> m_RdLookUpTable, m_TdLookUpTable;
};

MTS_IMPLEMENT_CLASS_S(IsotropicMultipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicMultipole, "Isotropic multipole model");
MTS_NAMESPACE_END
