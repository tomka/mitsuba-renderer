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

MTS_NAMESPACE_BEGIN

typedef SubsurfaceMaterialManager::LUTType LUTType;
typedef SubsurfaceMaterialManager::LUTRecord LUTRecord;

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources
 */
struct IsotropicDipoleQuery {
#if !defined(MTS_SSE) || (SPECTRUM_SAMPLES != 3)
	inline IsotropicDipoleQuery(const Spectrum &zr, const Spectrum &zv, 
		const Spectrum &sigmaTr, Float Fdt, const Point &p)
		: zr(zr), zv(zv), sigmaTr(sigmaTr), result(0.0f), Fdt(Fdt), p(p) {
			count = 0;
	}

	inline void operator()(const IrradianceSample &sample) {
		Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());
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
	}

	inline void operator()(const IrradianceSample &sample) {
		/* Distance to the positive point source of the dipole */
		const __m128 lengthSquared = _mm_set1_ps((p - sample.p).lengthSquared()),
			drSqr = _mm_add_ps(zrSqr, lengthSquared), 
			dvSqr = _mm_add_ps(zvSqr, lengthSquared),
			dr = _mm_sqrt_ps(drSqr), dv = _mm_sqrt_ps(dvSqr), 
			one = _mm_set1_ps(1.0f),
			factor = _mm_mul_ps(_mm_set1_ps(0.25f*INV_PI*sample.area * Fdt),
				_mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0)),
			C1fac = _mm_div_ps(_mm_mul_ps(zr, _mm_add_ps(sigmaTr, _mm_div_ps(one, dr))), drSqr),
			C2fac = _mm_div_ps(_mm_mul_ps(zr, _mm_add_ps(sigmaTr, _mm_div_ps(one, dv))), dvSqr);
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
	Float Fdt;
	Point p;
};

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources
 */
struct IsotropicMultipoleQuery {
//#if !defined(MTS_SSE) || (SPECTRUM_SAMPLES != 3)
	inline IsotropicMultipoleQuery(const std::vector<Spectrum> &zr, const std::vector<Spectrum> &zv, 
		const Spectrum &sigmaTr, Float Fdt, const Point &p, const Normal &_ns,
        const Float _d, const Float numExtraDipoles)
		: zrList(zr), zvList(zv), sigmaTr(sigmaTr), result(0.0f), Fdt(Fdt), p(p),
          ns(_ns), d(_d), numExtraDipoles(numExtraDipoles), one(1.0f), count(0) {
	}

	inline void operator()(const IrradianceSample &sample) {
        const Spectrum negSigmaTr = sigmaTr * (-1.0f);
		const Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());

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

//#else
//	inline IsotropicMUltipoleQuery(const Spectrum &_zr, const Spectrum &_zv, 
//		const Spectrum &_sigmaTr, Float Fdt, const Point &p, const Normal &_ns,
//        const Float _d) : Fdt(Fdt), p(p), ns(_ns), d(_d) {
//		zr = _mm_set_ps(_zr[0], _zr[1], _zr[2], 0);
//		zv = _mm_set_ps(_zv[0], _zv[1], _zv[2], 0);
//		sigmaTr = _mm_set_ps(_sigmaTr[0], _sigmaTr[1], _sigmaTr[2], 0);
//		zrSqr = _mm_mul_ps(zr, zr);
//		zvSqr = _mm_mul_ps(zv, zv);
//		result.ps = _mm_setzero_ps();
//		count = 0;
//	}
//
//	inline void operator()(const IrradianceSample &sample) {
//		/* Distance to the positive point source of the dipole */
//		const __m128 lengthSquared = _mm_set1_ps((p - sample.p).lengthSquared()),
//			drSqr = _mm_add_ps(zrSqr, lengthSquared), 
//			dvSqr = _mm_add_ps(zvSqr, lengthSquared),
//			dr = _mm_sqrt_ps(drSqr), dv = _mm_sqrt_ps(dvSqr), 
//			one = _mm_set1_ps(1.0f),
//			factor = _mm_mul_ps(_mm_set1_ps(0.25f*INV_PI*sample.area * Fdt),
//				_mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0)),
//			C1fac = _mm_div_ps(_mm_mul_ps(distR, _mm_add_ps(sigmaTr, _mm_div_ps(one, dr))), drSqr),
//			C2fac = _mm_div_ps(_mm_mul_ps(distV, _mm_add_ps(sigmaTr, _mm_div_ps(one, dv))), dvSqr);
//		SSEVector temp1(_mm_mul_ps(dr, sigmaTr)), temp2(_mm_mul_ps(dv, sigmaTr));
//		const __m128
//			exp1 = _mm_set_ps(expf(-temp1.f[3]), expf(-temp1.f[2]), expf(-temp1.f[1]), 0),
//			exp2 = _mm_set_ps(expf(-temp2.f[3]), expf(-temp2.f[2]), expf(-temp2.f[1]), 0);
//		result.ps = _mm_add_ps(result.ps, _mm_mul_ps(factor, _mm_add_ps(
//			_mm_mul_ps(C1fac, exp1), _mm_mul_ps(C2fac, exp2))));
//	}
//
//	Spectrum getResult() {
//		Spectrum value;
//		for (int i=0; i<3; ++i)
//			value[i] = result.f[3-i];
//		return value;
//	}
//
//	__m128 zr, zv, zrSqr, zvSqr, sigmaTr;
//	SSEVector result;
//#endif

	Float Fdt;
	Point p;
    Normal ns;
    Spectrum d;
    int numExtraDipoles;
    const Spectrum one;
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
            Float _res, Float _Fdt, const Point &_p, const Normal &_ns)
        : dMoR_LUT(lutR), dMoT_LUT(lutT), entries(lutR->size()), invResolution(1.0f / _res),
          result(0.0f), Fdt(_Fdt), p(_p), ns(_ns), count(0) {
	}

	inline void operator()(const IrradianceSample &sample) {
	    //Float rSqr = (p - sample.p).lengthSquared();
	    Float r = (p - sample.p).length();
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
};


static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme with his and Donner's multipole approach.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *   Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 *
 * ("Light Diffusion in Multi-Layered Tranclucent Materials", in ACM
 *   Transactions on Graphics 24)
 */
class IsotropicMultipole : public Subsurface {
public:
	IsotropicMultipole(const Properties &props) 
		: Subsurface(props) {
		irrOctreeMutex->lock();
		m_octreeIndex = irrOctreeIndex++;
		irrOctreeMutex->unlock();
		
		/* How many samples should be taken when estimating the irradiance at a given point in the scene? 
		This attribute is currently only used in conjunction with subsurface integrators and
		can safely be ignored if the scene contains none of them. */
		m_irrSamples = props.getInteger("irrSamples", 32);
				
		/* When estimating the irradiance at a given point, should indirect illumination be included
		in the final estimate? This attribute is currently only used in conjunction with 
		subsurface integrators and can safely be ignored if the scene contains none of them. */
		m_irrIndirect = props.getBoolean("irrIndirect", true);

		/* Multiplicative factor, which can be used to adjust the number of
		   irradiance samples */
		m_sampleMultiplier = props.getFloat("sampleMultiplier", 2.0f);
		/* Error threshold - lower means better quality */
		m_minDelta= props.getFloat("quality", 0.1f);
		/* Max. depth of the created octree */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Multiplicative factor for the subsurface term - can be used to remove
		   this contribution completely, making it possible to use this integrator
		   for other interesting things.. */
		m_ssFactor = props.getSpectrum("ssFactor", Spectrum(1.0f));
		m_maxDepth = props.getInteger("maxDepth", 40);
		m_extraDipoles = props.getInteger("extraDipoles", 0);
        m_slabThickness = props.getFloat("slabThickness", 1.0f);
        m_useMartelliD = props.getBoolean("useMartelliDC", true);
        /* look-up-table */
        m_useRdLookUpTable = props.getBoolean("useLookUpTable", true);
        m_errThreshold = props.getFloat("errorThreshold", 0.01);
        m_lutResolution = props.getFloat("lutResolution", 0.01);
        m_rMaxPredefined = props.hasProperty("lutRmax");
        if (m_rMaxPredefined) {
            if (props.hasProperty("mcIterations"))
                Log(EError, "You can either specify 'lutRMax' or 'mcIterations', not both.");

            m_rMax = props.getFloat("lutRmax");
        }
        m_mcIterations = props.getInteger("mcIterations", 10000);

        if (m_extraDipoles == 0) {
            Log(EInfo, "Using standard dipole model");
        } else {
            Log(EInfo, "Using multipole model with %d dipoles in total", 2 * m_extraDipoles + 1);
        }

		/* Asymmetry parameter of the phase function */
		m_g = props.getFloat("g", 0);
		m_ready = false;
		m_octreeResID = -1;
	}
	
	IsotropicMultipole(Stream *stream, InstanceManager *manager) 
	 : Subsurface(stream, manager) {
		m_ssFactor = Spectrum(stream);
		m_g = stream->readFloat();
		m_sampleMultiplier = stream->readFloat();
		m_minDelta = stream->readFloat();
		m_maxDepth = stream->readInt();
		m_octreeIndex = stream->readInt();
		m_irrSamples = stream->readInt();
		m_irrIndirect = stream->readBool();
        m_useMartelliD = stream->readBool();
        m_useRdLookUpTable = stream->readBool();
        m_errThreshold = stream->readFloat();
        m_lutResolution = stream->readFloat();
		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~IsotropicMultipole() {
		if (m_octreeResID != -1)
			Scheduler::getInstance()->unregisterResource(m_octreeResID);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		if (m_octreeResID != -1)
			proc->bindResource(formatString("irrOctree%i", m_octreeIndex), m_octreeResID);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		m_ssFactor.serialize(stream);
		stream->writeFloat(m_g);
		stream->writeFloat(m_sampleMultiplier);
		stream->writeFloat(m_minDelta);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
		stream->writeInt(m_octreeIndex);
        stream->writeBool(m_useMartelliD);
        stream->writeBool(m_useRdLookUpTable);
        stream->writeFloat(m_errThreshold);
        stream->writeFloat(m_lutResolution);
	}

	Spectrum Lo(const Scene *scene, const Intersection &its, const Vector &d) const {
		if (!m_ready || m_ssFactor.isZero())
			return Spectrum(0.0f);

        const Normal &n = its.geoFrame.n;

        Spectrum Mo = Spectrum(0.0f);
        /* If there are no extra dipoles, do a srandard dipolo query */
        if (m_extraDipoles == 0) {
                IsotropicDipoleQuery query(m_zr[0], m_zv[0], m_sigmaTr, m_Fdt, its.p);
                m_octree->execute(query);
                // compute combined radiant exitance
                Mo = query.getResult();
        } else {
            if (m_useRdLookUpTable) {
                IsotropicLUTMultipoleQuery query(m_RdLookUpTable, m_TdLookUpTable, m_lutResolution, m_Fdt, its.p, n);
                m_octree->execute(query);
                Mo = query.getResult();
            } else {
                // calulate diffuse reflectance and transmittance
                IsotropicMultipoleQuery query(m_zr, m_zv, m_sigmaTr, m_Fdt, its.p, n, m_slabThickness, m_extraDipoles);
                m_octree->execute(query);

                // compute combined radiant exitance
                Mo += query.getResult();
            }
        }

		if (m_eta == 1.0f) {
			return Mo * m_ssFactor * INV_PI;
		} else {
			Float Ft = 1.0f - fresnel(absDot(n, d));
			return Mo * m_ssFactor * INV_PI * (Ft / m_Fdr);
		}
	}

	void configure() {
		m_sigmaSPrime = m_sigmaS * (1-m_g);
		m_sigmaTPrime = m_sigmaSPrime + m_sigmaA;

		/* Mean-free path (avg. distance traveled through the medium) */
		m_mfp = Spectrum(1.0f) / m_sigmaTPrime;

		/* Also find the smallest mean-free path for all wavelengths */
		m_minMFP = std::numeric_limits<Float>::max();
		for (int lambda=0; lambda<SPECTRUM_SAMPLES; lambda++)
			m_minMFP = std::min(m_minMFP, m_mfp[lambda]);

        if (m_eta > 1) {
            /* Average reflectance due to mismatched indices of refraction
               at the boundary - [Groenhuis et al. 1983]*/
            m_Fdr = -1.440f / (m_eta * m_eta) + 0.710f / m_eta
                + 0.668f + 0.0636f * m_eta;
        } else {
            /* Average reflectance due to mismatched indices of refraction
               at the boundary - [Egan et al. 1973]*/
            m_Fdr = -0.4399f + 0.7099f / m_eta - 0.3319f / (m_eta * m_eta)
                + 0.0636f / (m_eta * m_eta * m_eta);
        }

		/* Reduced albedo */
		m_alphaPrime = m_sigmaSPrime / m_sigmaTPrime;

		/* Average transmittance at the boundary */
		m_Fdt = 1.0f - m_Fdr;

		if (m_eta == 1.0f) {
			m_Fdr = (Float) 0.0f;
			m_Fdt = (Float) 1.0f;
		}

		/* Approximate dipole boundary condition term */
		m_A = (1 + m_Fdr) / m_Fdt;

		/* Effective transport extinction coefficient */
		m_sigmaTr = (m_sigmaA * m_sigmaTPrime * 3.0f).sqrt();

		/* Diffusion coefficient
         * According to Martelli et al.'s paper "Accuracy of the Diffusion
         * Equation to Describe Photon Migration through an Infininite
         * Medium" from 2000, the diffusion coefficient should be calculated
         * slightly different. In practice this seems only required when
         * sigmaA / sigmaSPrime > 0.01.
         */
        if (m_useMartelliD)
		    m_D = Spectrum(1.0f) / (m_sigmaSPrime * 3.0f + m_sigmaA);
        else
		    m_D = Spectrum(1.0f) / (m_sigmaTPrime * 3.0f);

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

        /* Create look-up-table for dMoR (if requested) that is indexed
         * by the squared distance to the sample. The Monte Carlo
         * is currently done with 10k samples. This could be done more
         * dynamic. */
        if (m_useRdLookUpTable) {
            ref<SubsurfaceMaterialManager> smm = SubsurfaceMaterialManager::getInstance();
            std::string lutHashR = smm->getMultipoleLUTHashR(m_lutResolution, m_errThreshold,
                m_sigmaTr, m_alphaPrime, m_extraDipoles, m_zr, m_zv);
            std::string lutHashT = smm->getMultipoleLUTHashT(m_lutResolution, m_errThreshold,
                m_sigmaTr, m_alphaPrime, m_extraDipoles, m_zr, m_zv, m_slabThickness);
            if (smm->hasLUT(lutHashR) && smm->hasLUT(lutHashT)) {
                /* get Rd LUT */
                LUTRecord lutR = smm->getLUT(lutHashR);
                m_RdLookUpTable = lutR.lut;
                AssertEx(lutR.resolution == m_lutResolution, "Cached Rd LUT does not have requested resolution!");
                /* get Td LUT */
                LUTRecord lutT = smm->getLUT(lutHashT);
                m_TdLookUpTable = lutT.lut;
                AssertEx(lutT.resolution == m_lutResolution, "Cached Td LUT does not have requested resolution!");
            } else {
                boost::timer timer;
                if (!m_rMaxPredefined) {
                    const Spectrum invSigmaTr = 1.0f / m_sigmaTr;
                    const Float inv4Pi = 1.0f / (4 * M_PI);
                    ref<Random> random = new Random();

                    /* Find Rd for the whole area by monte carlo integration. The
                     * sampling area is calculated from the max. mean free path.
                     * A square area around with edge length 2 * maxMFP is used
                     * for this. Hence, the sampling area is 4 * maxMFP * maxMFP. */
                    const int numSamples = m_mcIterations;
                    Spectrum Rd_A = Spectrum(0.0f);
                    for (int n = 0; n < numSamples; ++n) {
                        /* do importance sampling by choosing samples distributed
                         * with sigmaTr^2 * e^(-sigmaTr * r). */
                        Spectrum r = invSigmaTr * -std::log( random->nextFloat() );
                        Rd_A += getRd(r);
                    }
                    Float A = 4 * invSigmaTr.max() * invSigmaTr.max();
                    Rd_A = A * Rd_A * m_alphaPrime * inv4Pi / (Float)(numSamples - 1);
                    Log(EDebug, "After %i MC integration iterations, Rd seems to be %s (took %.2fs)",
                        numSamples, Rd_A.toString().c_str(), timer.elapsed());

                    /* Since we now have Rd integrated over the whole surface we can find a valid rmax
                     * for the given threshold. */
                    timer.restart();
                    const Float step = m_lutResolution;
                    Float rMax = 0.0f;
                    Spectrum err(std::numeric_limits<Float>::max());
                    while (err.max() > m_errThreshold) {
                        rMax += step;
                        /* Again, do MC integration, but with r clamped at rmax. */
                        Spectrum Rd_APrime(0.0f);
                        for (int n = 0; n < numSamples; ++n) {
                            /* do importance sampling by choosing samples distributed
                             * with sigmaTr^2 * e^(-sigmaTr * r). */
                            Spectrum r = invSigmaTr * -std::log( random->nextFloat() );
                            // clamp samples to rMax
                            for (int s=0; s<SPECTRUM_SAMPLES; ++s) {
                                r[s] = std::min(rMax, r[s]);
                            }
                            Rd_APrime += getRd(r);
                        }
                        Float APrime = 4 * rMax * rMax;
                        Rd_APrime = APrime * Rd_APrime * m_alphaPrime * inv4Pi / (Float)(numSamples - 1);
                        err = (Rd_A - Rd_APrime) / Rd_A;
                    }
                    m_rMax = rMax;
                    Log(EDebug, "Maximum distance for sampling surface is %f with an error of %f (took %.0fs)",
                        m_rMax, m_errThreshold, timer.elapsed());
                }

                /* Create the actual look-up-table if it was MC integrated */
                timer.restart();
                const int numEntries = (int) (m_rMax / m_lutResolution) + 1;
                m_RdLookUpTable = new LUTType(numEntries);
                m_TdLookUpTable = new LUTType(numEntries);
                for (int i=0; i<numEntries; ++i) {
                    const Spectrum r(i * m_lutResolution);
                    m_RdLookUpTable->at(i) = getRd(r);
                    m_TdLookUpTable->at(i) = getTd(r);
                }

                /* Create new LUTRecord and store this LUT */
                LUTRecord lutRecR(m_RdLookUpTable, m_lutResolution);
                smm->addLUT(lutHashR, lutRecR);
                LUTRecord lutRecT(m_RdLookUpTable, m_lutResolution);
                smm->addLUT(lutHashT, lutRecT);
                AssertEx(smm->hasLUT(lutHashR), "Rd LUT is not available, but it should be!");
                AssertEx(smm->hasLUT(lutHashT), "Td LUT is not available, but it should be!");

                Log(EDebug, "Created Rd and Td look-up-table with %i entries each (took %.2fs)", numEntries, timer.elapsed());
            }
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

	/// Unpolarized fresnel reflection term for dielectric materials
	Float fresnel(Float cosThetaI) const {
		Float g = std::sqrt(m_eta*m_eta - 1.0f + cosThetaI * cosThetaI);
		Float temp1 = (g - cosThetaI)/(g + cosThetaI);
		Float temp2 = (cosThetaI * (g + cosThetaI) - 1) / 
			(cosThetaI * (g - cosThetaI) + 1.0f);
		return 0.5f * temp1 * temp1 * (1.0f + temp2 * temp2);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
		if (m_ready)
			return true;

		if (!scene->getIntegrator()->getClass()
				->derivesFrom(SampleIntegrator::m_theClass)) {
			Log(EError, "The dipole subsurface integrator requires "
				"a sampling-based surface integrator!");
		}

		m_octree = new IrradianceOctree(m_maxDepth, m_minDelta, 
			scene->getKDTree()->getAABB());

		Float sa = 0;
		for (std::vector<Shape *>::iterator it = m_shapes.begin(); 
			it != m_shapes.end(); ++it) 
			sa += (*it)->getSurfaceArea();
		size_t sampleCount = (size_t) std::ceil(sa / (M_PI * m_minMFP * m_minMFP)
			* m_sampleMultiplier);
		Log(EInfo, "Generating " SIZE_T_FMT " irradiance samples..", sampleCount);

		ref<Scheduler> sched = Scheduler::getInstance();

		/* This could be a bit more elegant.. - inform the irradiance
		   sampler about the index of this subsurface integrator */
		std::vector<Subsurface *> ssIntegrators
			= scene->getSubsurfaceIntegrators();
		int index = -1;
		for (size_t i=0; i<ssIntegrators.size(); ++i) {
			if (ssIntegrators[i] == this) {
				index = i;
				break;
			}
		}
		Assert(index != -1);

		ref<IrradianceSamplingProcess> proc = new IrradianceSamplingProcess(
			sampleCount, (size_t) std::ceil(sampleCount/100.0f), index, 
			m_irrSamples, m_irrIndirect, job);

		proc->bindResource("scene", sceneResID);
		scene->bindUsedResources(proc);
		m_proc = proc;
		sched->schedule(proc);
		sched->wait(proc);
		m_proc = NULL;
		if (proc->getReturnStatus() != ParallelProcess::ESuccess)
			return false;

		const IrradianceRecordVector &results = *proc->getSamples();
		for (size_t i=0; i<results.size(); ++i) 
			m_octree->addSample(results[i]);

		m_octree->preprocess();
		m_octreeResID = Scheduler::getInstance()->registerResource(m_octree);

		m_ready = true;
		return true;
	}

	void wakeup(std::map<std::string, SerializableObject *> &params) {
		std::string octreeName = formatString("irrOctree%i", m_octreeIndex);
		if (!m_octree.get() && params.find(octreeName) != params.end()) {
			m_octree = static_cast<IrradianceOctree *>(params[octreeName]);
			m_ready = true;
		}
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_proc);
	}

	MTS_DECLARE_CLASS()
private:
	Float m_minMFP, m_sampleMultiplier;
	Float m_Fdr, m_Fdt, m_A, m_minDelta, m_g;
	Spectrum m_mfp, m_sigmaTr, m_alphaPrime;
	Spectrum m_sigmaSPrime, m_sigmaTPrime, m_D, m_ssFactor;
    /* True if alternate diffusion constant calculation should be used */
    bool m_useMartelliD;
    /* Lists for the dipole distances */
    std::vector<Spectrum> m_zr, m_zv;
	ref<IrradianceOctree> m_octree;
	ref<ParallelProcess> m_proc;
	int m_octreeResID, m_octreeIndex;
	int m_maxDepth;
	int m_irrSamples;
	bool m_irrIndirect;
	bool m_ready, m_requireSample;
    /* The number of additional dipoles. Specifies the number of
     * dipoles that get added below *and* above the surface. E.g.
     * 1 additional dipole means you have the standard dipole plus
     * one above and one below.
     */
    int m_extraDipoles;
    /* The slab thickness */
    Float m_slabThickness;
    /* Indicates if a look-up-table should be created and used for Rd */
    bool m_useRdLookUpTable;
    /* Look-up-table for Rd and Td, indexed by the distance r. */
    ref<LUTType> m_RdLookUpTable, m_TdLookUpTable;
    /* the maximum distance stored in the LUT */
    Float m_rMax;
    /* is rMax predefined? */
    bool m_rMaxPredefined;
    /* error threshold for rmax */
    Float m_errThreshold;
    /* monte carlo integration iterations */
    int m_mcIterations;
    /* resolution of the dMoR LUT */
    Float m_lutResolution;
};

MTS_IMPLEMENT_CLASS_S(IsotropicMultipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicMultipole, "Isotropic multipole model");
MTS_NAMESPACE_END
