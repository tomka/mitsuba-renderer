/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2011 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <boost/bind.hpp>
#include <boost/timer.hpp>
#include "../medium/materials.h"
#include "irrtree.h"
#include "dipolebased.h"

//#define NO_SSE_QUERY

MTS_NAMESPACE_BEGIN

/* Relative bound on what is still accepted as roundoff 
   error -- be quite tolerant */
#if defined(SINGLE_PRECISION)
	#define ERROR_REQ 1e-2f
#else
	#define ERROR_REQ 1e-5
#endif

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources
 */
struct IsotropicDipoleQuery {
#if !defined(MTS_SSE) || (SPECTRUM_SAMPLES != 3) || defined(NO_SSE_QUERY)
	inline IsotropicDipoleQuery(const Spectrum &zr, const Spectrum &zv, 
		const Spectrum &sigmaTr, Float Fdt, const Point &p) 
		: zr(zr), zv(zv), sigmaTr(sigmaTr), result(0.0f), Fdt(Fdt), p(p) {
			count = 0;
			Float zrMin = zr.min();
			zrMinSq = zrMin * zrMin;
	}

	inline void operator()(const IrradianceSample &sample) {
		Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq); 
		//Float dist = (p - sample.p).lengthSquared(); 
		Spectrum rSqr = Spectrum(dist);
		/* Distance to the real source */
		Spectrum dr = (rSqr + zr*zr).sqrt();
		/* Distance to the image point source */
		Spectrum dv = (rSqr + zv*zv).sqrt();
		Spectrum C1 = zr * (sigmaTr + Spectrum(1.0f) / dr);
		Spectrum C2 = zv * (sigmaTr + Spectrum(1.0f) / dv);

		/* Do not include the reduced albedo - will be canceled out later */
		Spectrum dMo = Spectrum(0.25f * INV_PI) *
			 (C1 * ((-sigmaTr * dr).exp()) / (dr * dr)
			+ C2 * ((-sigmaTr * dv).exp()) / (dv * dv));
		result += dMo * sample.E * (sample.area * Fdt);
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	Spectrum zr, zv, sigmaTr, result;
#else
	inline IsotropicDipoleQuery(const Spectrum &_zr, const Spectrum &_zv, 
		const Spectrum &_sigmaTr, Float Fdt, const Point &p) : Fdt(Fdt), p(p) {
		zr = _mm_set_ps(_zr[0], _zr[1], _zr[2], 0);
		zv = _mm_set_ps(_zv[0], _zv[1], _zv[2], 0);
		sigmaTr = _mm_set_ps(_sigmaTr[0], _sigmaTr[1], _sigmaTr[2], 0);
		zrSqr = _mm_mul_ps(zr, zr);
		zvSqr = _mm_mul_ps(zv, zv);
		result.ps = _mm_setzero_ps();
		count = 0;
		Float zrMin = _zr.min();
		zrMinSq = zrMin * zrMin;
	}

	inline void operator()(const IrradianceSample &sample) {
		/* Distance to the positive point source of the dipole */
		Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq);
		const __m128 lengthSquared = _mm_set1_ps(dist),
			drSqr = _mm_add_ps(zrSqr, lengthSquared), 
			dvSqr = _mm_add_ps(zvSqr, lengthSquared),
			dr = _mm_sqrt_ps(drSqr), dv = _mm_sqrt_ps(dvSqr), 
			one = _mm_set1_ps(1.0f),
			factor = _mm_mul_ps(_mm_set1_ps(0.25f*INV_PI*sample.area * Fdt),
				_mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0)),
			C1fac = _mm_div_ps(_mm_mul_ps(zr, _mm_add_ps(sigmaTr, _mm_div_ps(one, dr))), drSqr),
			C2fac = _mm_div_ps(_mm_mul_ps(zv, _mm_add_ps(sigmaTr, _mm_div_ps(one, dv))), dvSqr);
		SSEVector temp1(_mm_mul_ps(dr, sigmaTr)), temp2(_mm_mul_ps(dv, sigmaTr));
		const __m128
			exp1 = _mm_set_ps(expf(-temp1.f[3]), expf(-temp1.f[2]), expf(-temp1.f[1]), 0),
			exp2 = _mm_set_ps(expf(-temp2.f[3]), expf(-temp2.f[2]), expf(-temp2.f[1]), 0);
		result.ps = _mm_add_ps(result.ps, _mm_mul_ps(factor, _mm_add_ps(
			_mm_mul_ps(C1fac, exp1), _mm_mul_ps(C2fac, exp2))));
	}

	Spectrum getResult() {
		Spectrum value;
		for (int i=0; i<3; ++i)
			value[i] = result.f[3-i];
		return value;
	}

	__m128 zr, zv, zrSqr, zvSqr, sigmaTr;
	SSEVector result;
#endif

	int count;
	Float Fdt, zrMinSq;
	Point p;
};

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources. This variant
 * requires a look-up-table and does currently not
 * benefit from SSE2.
 */
struct IsotropicLUTDipoleQuery {
#if !defined(MTS_SSE) || (SPECTRUM_SAMPLES != 3) || defined(NO_SSE_QUERY)
	inline IsotropicLUTDipoleQuery(const ref<LUTType> &lut, Float _res,
			Float Fdt, const Point &p, Float _minDist)
		: result(0.0f), dMo_LUT(lut), entries(lut->size()),
		  invResolution(1.0f / _res), Fdt(Fdt), p(p), count(0),
		  minDist(_minDist) { }

	inline void operator()(const IrradianceSample &sample) {
		Float r = (p - sample.p).length();
		// Avoid singularities (see Jensen et al. 2001)
		r = std::max(r, minDist);
		/* Look up dMo for the distance. As in the normal query,the
		 * reduced albedo is not included. It will be canceled out
		 * later. The index is rounded to the next nearest integer. */
		int index = (int) (r * invResolution + 0.5f);
		if (index < entries) {
			Spectrum dMo = dMo_LUT->at(index);
			result += dMo * sample.E * (sample.area * Fdt);
		}
		//TODO: Find out why we don't multiply with 0.25 * INV_PI
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	Spectrum result;
#else
	inline IsotropicLUTDipoleQuery(const ref<LUTType> &lut, Float _res,
			Float Fdt, const Point &p, Float _minDist)
		: dMo_LUT(lut), entries(lut->size()), invResolution(1.0f / _res),
		  Fdt(Fdt), p(p), count(0), minDist(_minDist) {
		result.ps = _mm_setzero_ps();
	}

	inline void operator()(const IrradianceSample &sample) {
		/* Distance to the positive point source of the dipole. Try
		 *Avoid singularities (see Jensen et al. 2001)
		 */
		Float dist = std::max((p - sample.p).length(), minDist);
		/* Look up dMo for the distance. As in the normal query,the
		 * reduced albedo is not included. It will be canceled out
		 * later. The index is rounded to the next nearest integer.
		 */
		int index = (int) (dist * invResolution + 0.5f);
		if (index < entries) {
			Spectrum dMo = dMo_LUT->at(index);

			const __m128 factor = _mm_mul_ps(_mm_set1_ps(sample.area * Fdt),
				_mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0));
			result.ps = _mm_add_ps(result.ps, _mm_mul_ps(factor,
							 _mm_set_ps(dMo[0], dMo[1], dMo[2], 0)));
		}
		count++;
	}

	Spectrum getResult() {
		Spectrum value;
		for (int i=0; i<3; ++i)
			value[i] = result.f[3-i];
		return value;
	}

	SSEVector result;
#endif
	/* LUT related */
	const ref<LUTType> &dMo_LUT;
	int entries;
	Float invResolution;
	
	//Float zrMinSq;
	Float Fdt;
	Point p;
	int count;
	Float minDist;
};

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *	 Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 */
class IsotropicDipole : public DipoleBasedSubsurface<IsotropicDipole> {
protected:

public:
	IsotropicDipole(const Properties &props) 
		: DipoleBasedSubsurface(props) { }

	IsotropicDipole(Stream *stream, InstanceManager *manager) 
		: DipoleBasedSubsurface(stream, manager) { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		DipoleBasedSubsurface::serialize(stream, manager);
	}

	Spectrum Mo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		IsotropicDipoleQuery query(m_zr, m_zv, m_sigmaTr, m_Fdt, its.p);
		m_octree->execute(query);

		return query.getResult();
	}

	Spectrum MoLUT(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		IsotropicLUTDipoleQuery query(m_RdLookUpTable, m_lutResolution, m_Fdt, its.p, m_minMFP);
		m_octree->execute(query);

		return query.getResult();
	}

	void configureDipoles() {
		/* Distance of the dipole point sources to the surface */
		m_zr = m_mfp; 
		m_zv = m_mfp * (1.0f + 4.0f/3.0f * m_A);
	}

	/* Creates a hash for the current configuraiton of this isotropic
	 * dipole. All relevant settings are included. This could e.g. be
	 * used to identify a cached look-up-table.
	 */
	std::string hashConfig() {
		std::ostringstream oss;
		// set precision to get around floating point rounding issues
		oss.precision(5);
		oss << "multipole" << m_lutResolution << "," << m_errThreshold << "," << m_sigmaTr.toString() << ","
			<< m_alphaPrime.toString() << "," << m_zr.toString() << "," << m_zv.toString();

		return oss.str();
	}

	void configureLUT() {
		ref<SubsurfaceMaterialManager> smm =
			 SubsurfaceMaterialManager::getInstance();
		std::string lutHash = hashConfig();

		if (smm->hasLUT(lutHash)) {
			LUTRecord lutR = smm->getLUT(lutHash);
			m_RdLookUpTable = lutR.lut;
			AssertEx(lutR.resolution == m_lutResolution,
				"Cached LUT does not have requested resolution!");
		} else {
			// call createLUT() of base
			m_RdLookUpTable = createLUT(boost::bind(
				&mitsuba::IsotropicDipole::getRd, this, _1));

			/* Create new LUTRecord and store this LUT */
			LUTRecord lutRecR(m_RdLookUpTable, m_lutResolution);
			smm->addLUT(lutHash, lutRecR);
			/* Make sure the LUT got actually stored. */
			AssertEx(smm->hasLUT(lutHash), "Rd LUT is not available, but it should be!");
		}
	}

	/* Calculate Rd based on all dipoles and the requested distance
	 */
	Spectrum getRd(const Spectrum r) {
		//Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq); 
		const Spectrum rSqr = r * r;
		const Spectrum one(1.0f);
		const Spectrum negSigmaTr = m_sigmaTr * (-1.0f);

		/* Distance to the real source */
		Spectrum dr = (rSqr + m_zr*m_zr).sqrt();
		/* Distance to the image point source */
		Spectrum dv = (rSqr + m_zv*m_zv).sqrt();

		Spectrum C1 = m_zr * (m_sigmaTr + one / dr);
		Spectrum C2 = m_zv * (m_sigmaTr + one / dv);

		/* Do not include the reduced albedo - will be canceled out later */
		Spectrum dMo = Spectrum(0.25f * INV_PI) *
			 (C1 * ((negSigmaTr * dr).exp()) / (dr * dr)
			+ C2 * ((negSigmaTr * dv).exp()) / (dv * dv));

		dMo.clampNegative();
		return dMo;
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_zr, m_zv;
	/* Look-up-table for Rd, indexed by the distance r. */
	ref<LUTType> m_RdLookUpTable;
};

MTS_IMPLEMENT_CLASS_S(IsotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicDipole, "Isotropic dipole model");
MTS_NAMESPACE_END
