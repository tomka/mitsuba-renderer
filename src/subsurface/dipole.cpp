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
            Float zrMin = zr.min();
            zrMinSq = zrMin * zrMin;
	}

	inline void operator()(const IrradianceSample &sample) {
        Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq); 
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
    Float zrMinSq;
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
	inline IsotropicLUTDipoleQuery(const ref<LUTType> &lut, Float _res, 
		    Float Fdt, const Point &p) 
		: dMo_LUT(lut), entries(lut->size()), invResolution(1.0f / _res),
          result(0.0f), Fdt(Fdt), p(p), count(0) {
	}

	inline void operator()(const IrradianceSample &sample) {
        //Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq);
	    Float r = (p - sample.p).length();
        /* Look up dMo for the distance. As in the normal query,
         * the reduced albedo is not included. It will be canceled
         * out later. */
        int index = (int) (r * invResolution);
        if (index < entries) {
            Spectrum dMo = dMo_LUT->at(index);

            /* combine Mo based on R and Mo based on T to a new
             * Mo based on a combined profile P. */
            result += dMo * sample.E * (sample.area * Fdt);
 
    		count++;
        }
	}

	inline const Spectrum &getResult() const {
		return result;
	}

    /* LUT related */
    const ref<LUTType> &dMo_LUT;
    int entries;
    Float invResolution;
    
    //Float zrMinSq;
    Spectrum result;
	Float Fdt;
	Point p;
	int count;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *   Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 */
class IsotropicDipole : public Subsurface {
public:
	IsotropicDipole(const Properties &props) 
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
		m_minDelta = props.getFloat("quality", 0.1f);
		/* Max. depth of the created octree */
		m_maxDepth = props.getInteger("maxDepth", 40);
        /* Single scattering term */
        m_singleScattering = props.getBoolean("singleScattering", false);
        /* Should the irrtree be dumped? */
        m_dumpIrrtree = props.getBoolean("dumpIrrtree", false);
        m_dumpIrrtreePath = props.getString("dumpIrrtreePath", "");
		/* Multiplicative factor for the subsurface term - can be used to remove
		   this contribution completely, making it possible to use this integrator
		   for other interesting things.. */
		m_ssFactor = props.getSpectrum("ssFactor", Spectrum(1.0f));
		/* Asymmetry parameter of the phase function */
		m_g = props.getFloat("g", 0);
        /* alternative diffusion coefficient */
        m_useMartelliD = props.getBoolean("useMartelliDC", true);
        /* texture usage */
        m_useTextures = props.getBoolean("useTexture", false);
        if (m_useTextures) {
            fs::path filename = Thread::getThread()->getFileResolver()->resolve(
                    props.getString("zrFilename"));
            Log(EInfo, "Loading texture \"%s\"", filename.leaf().c_str());

            ref<FileStream> fs = new FileStream(filename, FileStream::EReadOnly);
            m_zrBitmap = new Bitmap(Bitmap::EEXR, fs);

            filename = Thread::getThread()->getFileResolver()->resolve(
                    props.getString("sigmaTrFilename"));
            Log(EInfo, "Loading texture \"%s\"", filename.leaf().c_str());

            fs = new FileStream(filename, FileStream::EReadOnly);
            m_sigmaTrBitmap = new Bitmap(Bitmap::EEXR, fs);

            m_texUScaling = props.getFloat("texUScaling", 1.0f);
            m_texVScaling = props.getFloat("texVScaling", 1.0f);
        }

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

		m_ready = false;
		m_octreeResID = -1;
	}
	
	IsotropicDipole(Stream *stream, InstanceManager *manager) 
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
        m_useTextures = stream->readBool();
        m_useRdLookUpTable = stream->readBool();
        m_errThreshold = stream->readFloat();
        m_lutResolution = stream->readFloat();
        m_mcIterations = stream->readInt();
		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~IsotropicDipole() {
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
		stream->writeInt(m_octreeIndex);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
        stream->writeBool(m_useMartelliD);
        stream->writeBool(m_useTextures);
        stream->writeBool(m_useRdLookUpTable);
        stream->writeFloat(m_errThreshold);
        stream->writeFloat(m_lutResolution);
        stream->writeInt(m_mcIterations);
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		if (!m_ready || m_ssFactor.isZero())
			return Spectrum(0.0f);

        if (m_useTextures) {
            Spectrum zr = m_zrTex->getValue(its);
            Spectrum zv = m_zvTex->getValue(its);
            Spectrum sigmaTr = m_sigmaTrTex->getValue(its);

            IsotropicDipoleQuery query(zr, zv, sigmaTr, m_Fdt, its.p);
            m_octree->execute(query);
            // compute multiple scattering term
            Spectrum Mo = query.getResult();

            const Normal &n = its.shFrame.n;
            Spectrum Lo;
            if (m_eta == 1.0f) {
                Lo = Mo * m_ssFactor * INV_PI;
            } else {
                Float Ft = 1.0f - fresnel(absDot(n, d));
                Lo = Mo * m_ssFactor * INV_PI * (Ft / m_Fdr);
            }
            return Lo;
        } else {
            Spectrum Mo;
            if (m_useRdLookUpTable) {
                IsotropicLUTDipoleQuery query(m_RdLookUpTable, m_lutResolution, m_Fdt, its.p);
                m_octree->execute(query);
                // compute multiple scattering term
                Mo = query.getResult();
            } else {
                IsotropicDipoleQuery query(m_zr, m_zv, m_sigmaTr, m_Fdt, its.p);
                m_octree->execute(query);
                // compute multiple scattering term
                Mo = query.getResult();
            }
        
            const Normal &n = its.shFrame.n;
            Spectrum Lo;
            if (m_eta == 1.0f) {
                Lo = Mo * m_ssFactor * INV_PI;
            } else {
                Float Ft = 1.0f - fresnel(absDot(n, d));
                Lo = Mo * m_ssFactor * INV_PI * (Ft / m_Fdr);
            }

            /* Compute single scattering term if requested. This is done with
             * one shadow ray per light. Then per shadow ray a number of samples
             * is used to calculate the contribution due to that one. */
            if (m_singleScattering) {
                const int nrSamples = 5;
                Float singleScatteringLo = 0.0f;
                for (int i=0; i < nrSamples; ++i) {
                    Vector wo;
                    //singleScatteringLo += LoSingleScattering(wo, d, its);
                }

                //Lo += singleScatteringLo / nrSamples;
            }

            return Lo;
        }
	}

	void configure() {
		m_sigmaSPrime = m_sigmaS * (1-m_g);
		m_sigmaTPrime = m_sigmaSPrime + m_sigmaA;

        /* extinction coefficient */
        m_sigmaT = m_sigmaA + m_sigmaS;
        /* get the longest mean free path */
        m_invSigmaTMin = 1.0f / m_sigmaT.min();
        m_invSigmaT = m_sigmaT.pow(-1.0f);
        m_negSigmaT = m_sigmaT * (-1.0f);

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

		/* Distance of the dipole point sources to the surface */
		m_zr = m_mfp; 
		m_zv = m_mfp * (1.0f + 4.0f/3.0f * m_A);

        /* Configure bitmap usage */
        if (m_useTextures)
            configureTexture();

        /* Configure look-up-table */
        if (m_useRdLookUpTable) {
            ref<SubsurfaceMaterialManager> smm = SubsurfaceMaterialManager::getInstance();
            std::string lutHash = smm->getDipoleLUTHash(m_lutResolution, m_errThreshold,
                m_sigmaTr, m_alphaPrime, m_zr, m_zv);
            if (smm->hasLUT(lutHash)) {
                LUTRecord lutR = smm->getLUT(lutHash);
                m_RdLookUpTable = lutR.lut;
                AssertEx(lutR.resolution == m_lutResolution, "Cached LUT does not have requested resolution!");
            } else {
                if (!m_rMaxPredefined) {
                    const Spectrum invSigmaTr = 1.0f / m_sigmaTr;
                    const Float inv4Pi = 1.0f / (4 * M_PI);
                    ref<Random> random = new Random();

                    /* Find Rd for the whole area by monte carlo integration. The
                     * sampling area is calculated from the max. mean free path.
                     * A square area around with edge length 2 * maxMFP is used
                     * for this. Hene, the sampling area is 4 * maxMFP * maxMFP. */
                    Spectrum Rd_A = Spectrum(0.0f);
                    int count;
                    for (count = 0; count < m_mcIterations; ++count) {
                        /* do importance sampling by choosing samples distributed
                         * with sigmaTr^2 * e^(-sigmaTr * r). */
                        Spectrum r = invSigmaTr * -std::log( random->nextFloat() );
                        Rd_A += getRd(r);
                    }
                    Float A = 4 * invSigmaTr.max() * invSigmaTr.max();
                    Rd_A = A * Rd_A * m_alphaPrime * inv4Pi / (Float)(m_mcIterations - 1);
                    Log(EDebug, "After %i MC integration iterations, Rd seems to be %s", count, Rd_A.toString().c_str());

                    /* Since we now have Rd integrated over the whole surface we can find a valid rmax
                     * for the given threshold. */
                    Float rMax = 0.0f;
                    Spectrum err(std::numeric_limits<Float>::max());
                    Spectrum invRd_A = Spectrum(1.0f) / Rd_A;
                    while (err.max() > m_errThreshold) {
                        rMax += m_lutResolution;
                        /* Again, do MC integration, but with r clamped at rmax. */
                        Spectrum Rd_APrime(0.0f);
                        for (int n = 0; n < m_mcIterations; ++n) {
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
                        Rd_APrime = APrime * Rd_APrime * m_alphaPrime * inv4Pi / (Float)(m_mcIterations - 1);
                        err = (Rd_A - Rd_APrime) * invRd_A;
                    }
                    m_rMax = rMax;
                    Log(EDebug, "Maximum distance for sampling surface is %f with an error of %f", m_rMax, m_errThreshold);
                }

                /* Create the actual look-up-table */
                const int numEntries = (int) (m_rMax / m_lutResolution) + 1;
                m_RdLookUpTable = new LUTType(numEntries);
                for (int i=0; i<numEntries; ++i) {
                    m_RdLookUpTable->at(i) = getdMoR(i * m_lutResolution);
                }

                /* Create new LUTRecord and store this LUT if it was MC integrated */
                if (!m_rMaxPredefined) {
                    LUTRecord lutRec(m_RdLookUpTable, m_lutResolution);
                    smm->addLUT(lutHash, lutRec);
                    AssertEx(smm->hasLUT(lutHash), "LUT is not available, but it should be!");
                }

                Log(EDebug, "Created Rd look-up-table with %i entries.", numEntries);
            }
        }
    }

    void configureTexture() {
        Random *random = new Random();
        m_random.set(random);

        PluginManager *pluginManager = PluginManager::getInstance();
        int w = m_zrBitmap->getWidth();
        int h = m_zrBitmap->getHeight();
        float *data = m_zrBitmap->getFloatData();

        /* create zr bitmap */
        ref<Bitmap> zrBitmap = new Bitmap(w, h, 128);
        float *zrData = zrBitmap->getFloatData();
       
        const bool adjustMFP = true;
        const Float origMinMFP = m_minMFP;

        /* if alpha of image is > 0, then use the RGB values */
        for (int y=0; y<m_zrBitmap->getHeight(); ++y) {
            for (int x=0; x<m_zrBitmap->getWidth(); ++x) {
                float r = *data++;
                float g = *data++;
                float b = *data++;
                ++data; // alpha
                float sum = r + g + b;

                if (sum > 0.001) {
                    *zrData++ = r;
                    *zrData++ = g;
                    *zrData++ = b;
                    *zrData++ = 1.0; // alpha
                    // find a potentially lower MFP
                    if (adjustMFP) {
                        /* the tests are rearranged for faster computation */
                        if (r < m_minMFP)
                            m_minMFP = r;
                        if (g < m_minMFP)
                            m_minMFP = g;
                        if (b < m_minMFP)
                            m_minMFP = b;
                    }
                } else {
                    *zrData++ = m_zr[0];
                    *zrData++ = m_zr[1];
                    *zrData++ = m_zr[2];
                    *zrData++ = 1.0;
                }
            }
        }
        /* write out the bitmap */
        std::string zrFileName = "zr" + randomString(m_random.get(), 7) + ".exr";
        fs::path filename = Thread::getThread()->getFileResolver()->resolve(zrFileName);

        Log(EInfo, "Writing zr texture \"%s\"", filename.leaf().c_str());

        ref<FileStream> outStream = new FileStream(filename, FileStream::ETruncWrite);
        zrBitmap->save(Bitmap::EEXR, outStream);
        outStream->close();

        /* create zr texture */
        Properties props;
        props.setPluginName("diffusiontexture");
        props.setString("filename", zrFileName);
        props.setFloat("uscale", m_texUScaling);
        props.setFloat("vscale", m_texVScaling);
        m_zrTex = static_cast<Texture *> (pluginManager->createObject(
                Texture::m_theClass, props));

        /* create zv bitmap */
        ref<Bitmap> zvBitmap = new Bitmap(w, h, 128);
        float *zvData = zvBitmap->getFloatData();
        zrData = zrBitmap->getFloatData();

        /* if alpha of image is > 0, then use the RGB values */
        for (int y=0; y<m_zrBitmap->getHeight(); ++y) {
            for (int x=0; x<m_zrBitmap->getWidth(); ++x) {
                float r = *zrData++;
                float g = *zrData++;
                float b = *zrData++;
                ++zrData; // alpha

                *zvData++ = r * (1.0f + (4.0f/3.0f) * m_A);
                *zvData++ = g * (1.0f + (4.0f/3.0f) * m_A);
                *zvData++ = b * (1.0f + (4.0f/3.0f) * m_A);
                *zvData++ = 1.0;
            }
        }
        /* write out the bitmap */
        std::string zvFileName = "zv" + randomString(m_random.get(), 7) + ".exr";
        filename = Thread::getThread()->getFileResolver()->resolve(zvFileName);

        Log(EInfo, "Writing zv texture \"%s\"", filename.leaf().c_str());

        outStream = new FileStream(filename, FileStream::ETruncWrite);
        zrBitmap->save(Bitmap::EEXR, outStream);
        outStream->close();

        /* create zv texture */
        Properties zvProps;
        zvProps.setPluginName("diffusiontexture");
        zvProps.setString("filename", zvFileName);
        zvProps.setFloat("uscale", m_texUScaling);
        zvProps.setFloat("vscale", m_texVScaling);
        m_zvTex = static_cast<Texture *> (pluginManager->createObject(
                Texture::m_theClass, zvProps));

        /* create sigmaTr bitmap */
        ref<Bitmap> sTrBitmap = new Bitmap(w, h, 128);
        float *sTrData = sTrBitmap->getFloatData();
        data = m_sigmaTrBitmap->getFloatData();

        /* if alpha of image is > 0, then use the RGB values */
        for (int y=0; y<m_sigmaTrBitmap->getHeight(); ++y) {
            for (int x=0; x<m_sigmaTrBitmap->getWidth(); ++x) {
                float r = *data++;
                float g = *data++;
                float b = *data++;
                ++data; // alpha
                float sum = r + g + b;

                if (sum > 0.001) {
                    *sTrData++ = r;
                    *sTrData++ = g;
                    *sTrData++ = b;
                    *sTrData++ = 1.0; // alpha
                } else {
                    *sTrData++ = m_sigmaTr[0];
                    *sTrData++ = m_sigmaTr[1];
                    *sTrData++ = m_sigmaTr[2];
                    *sTrData++ = 1.0;
                }
            }
        }
        /* write out the bitmap */
        std::string sTrFileName = "sigmaTr" + randomString(m_random.get(), 7) + ".exr";
        filename = Thread::getThread()->getFileResolver()->resolve(sTrFileName);

        Log(EInfo, "Writing sigmaTr texture \"%s\"", filename.leaf().c_str());

        outStream = new FileStream(filename, FileStream::ETruncWrite);
        sTrBitmap->save(Bitmap::EEXR, outStream);
        outStream->close();

        /* create sigmaTr texture */
        Properties sTrProps;
        sTrProps.setPluginName("diffusiontexture");
        sTrProps.setString("filename", sTrFileName);
        sTrProps.setFloat("uscale", m_texUScaling);
        sTrProps.setFloat("vscale", m_texVScaling);
        m_sigmaTrTex = static_cast<Texture *> (pluginManager->createObject(
                Texture::m_theClass, sTrProps));

        if (std::abs(origMinMFP - m_minMFP) > 0.0001) {
            Log(EInfo, "Adjusted minimum MFP from %.6f to %.6f", origMinMFP, m_minMFP);
        }
	}

    /// Calculate Rd based on all dipoles and the requested distance
    Spectrum getRd(Spectrum r) {
        const Spectrum one(1.0f);
        const Spectrum negSigmaTr = m_sigmaTr * (-1.0f);
		const Spectrum rSqr = r * r;

        // calulate diffuse reflectance and transmittance
        Spectrum dr = (rSqr + m_zr*m_zr).sqrt();
        Spectrum dv = (rSqr + m_zv*m_zv).sqrt();

        // the change in Rd
        Spectrum Rd =   (m_zr * (one + m_sigmaTr * dr) * (negSigmaTr * dr).exp() / (dr * dr * dr))
                      + (m_zv * (one + m_sigmaTr * dv) * (negSigmaTr * dv).exp() / (dv * dv * dv));
        return Rd;
    }

    Spectrum getdMoR(Float r) {
        //Float dist = std::max((p - sample.p).lengthSquared(), zrMinSq); 
		Spectrum rSqr = Spectrum(r * r);

        /* Distance to the real source */
        Spectrum dr = (rSqr + m_zr*m_zr).sqrt();
        /* Distance to the image point source */
        Spectrum dv = (rSqr + m_zv*m_zv).sqrt();

        Spectrum C1 = (m_sigmaTr + Spectrum(1.0f) / dr);
        Spectrum C2 = (m_sigmaTr + Spectrum(1.0f) / dv);

        /* Do not include the reduced albedo - will be canceled out later */
        Spectrum dMo = Spectrum(0.25f * INV_PI) *
             (m_zr * C1 * ((-m_sigmaTr * dr).exp()) / (dr * dr)
            + m_zv * C2 * ((-m_sigmaTr * dv).exp()) / (dv * dv));
        return dMo;
    }

	/// Unpolarized fresnel reflection term for dielectric materials
	Float fresnel(Float cosThetaI) const {
		Float g = std::sqrt(m_eta*m_eta - 1.0f + cosThetaI * cosThetaI);
		Float temp1 = (g - cosThetaI)/(g + cosThetaI);
		Float temp2 = (cosThetaI * (g + cosThetaI) - 1) / 
			(cosThetaI * (g - cosThetaI) + 1.0f);
		return 0.5f * temp1 * temp1 * (1.0f + temp2 * temp2);
	}

    /**
     * Computes the single-scattering radiance with the help of a BSSRDF.
     * ToDo: Actual Monte Carlo sampling.
     */
    Spectrum LoSingleScattering(const Vector &wi, const Vector &wo, const Intersection &its) const {
        /* cosines of input and output directions */
        const Float cos_wi = Frame::cosTheta(wi);
        const Float cos_wo = Frame::cosTheta(wo);
        const Float cos_wo_abs = std::abs(cos_wo);

        //Float eta = m_etaInt / m_etaExt;
        Float oneovereta = 1.0 / m_eta;
        Float oneoveretaSq = oneovereta * oneovereta;;

        /* Using Snell's law, calculate the squared sine of the
         * angle between the normal and the transmitted ray */
        Float sinTheta2Sqr = oneoveretaSq * Frame::sinTheta2(wi);

        if (sinTheta2Sqr > 1.0f) /* Total internal reflection! */
            return Spectrum(1.0f);

        /* Compute the cosine, but guard against numerical imprecision */
        Float cosTheta2 = std::sqrt(std::max((Float) 0.0f, 1.0f - sinTheta2Sqr));
        /* With cos(N, transmittedRay) on tap, calculating the 
         * transmission direction is straightforward. */
        Vector localTo = Vector(-oneovereta*wi.x, -oneovereta*wi.y, -cosTheta2);
        Vector to = normalize( its.toWorld( localTo ));

        /* importance sampling norminator */
        Random* random = m_random.get();
        if (random == NULL) {
            random = new Random();
            m_random.set(random);
        }
        Float sample = random->nextFloat();
        if (sample < 0.001)
            sample = 0.001;
        const Float ran = - std::log( sample );
        /* so' with max. maen frea path */
        const Float soPrimeMin =  m_invSigmaTMin * ran;

        /* Get sample point on refracted ray in world coordinates */
        const Point &xi = its.p;
        const Point3 xsamp = xi + to * soPrimeMin;

        /* Calculate siPrime and soPrime */

        /* Indireclty find intersection of light with surface xo by using the
         * triangle xi, xo, xamp with angles ai, ao, asamp. By using the
         * height/z-difference between xi and xsamp, we can calculate
         * si = h/(sin ao). si is the distance from sample point in surface
         * to light entering point. If gamma is the angle between normal an wo,
         * then ao = 90 degree - gamma. The sine of ao eqals the sine of
         * (90 - gamma), which again is (sin 90 * cos gamma - cos 90 * sin gamma)
         * This can be reduced to cos gamma.
         */
        const Float si = std::abs(xi.z - xsamp.z) / cos_wo;

        /* so' over whole spectrum */
        const Float term = 1.0f - (cos_wo_abs * cos_wo_abs);
        const Float siPrime = si * cos_wo_abs / sqrt(1.0f - oneoveretaSq * term);

        /* Calculate combined transmission coefficient */
        const Float G = std::abs(cosTheta2) / cos_wo_abs;
        const Spectrum sigmaTc = m_sigmaT + m_sigmaT * G;

        /* Calculate Fresnel trensmission T = 1- R */
        const Float Ft1 = 1.0f - mitsuba::fresnel(cos_wo, 1.0f, m_eta);
        const Float Ft2 = 1.0f - mitsuba::fresnel(cos_wi, 1.0f, m_eta);
        const Float F = Ft1 * Ft2;

        /* Query phase function */
        const Float p = hgPhaseFunction(wi, wo, m_g);

        const Spectrum siTerm = (m_negSigmaT * siPrime).exp();
        /* Actually the soTerm would be e^(-sPrime_o * sigmaT), but
         * this could be reduced to e^(ran) since sPrime_o = -ran/sigmaT. */
        const Spectrum soTerm = Spectrum( exp(ran) );

        Spectrum Lo = (m_sigmaS * F * p / sigmaTc) * siTerm * soTerm;
        return Lo;
    }

    /**
     *  Evaluate the Henyey-Greenstein phase function for two vectors with
        an asymmetry value g.  v1 and v2 should be normalized and g should 
     *  be in the range (-1, 1).  Negative values of g correspond to more
     *  back-scattering and positive values correspond to more forward scattering.
     */
    Float hgPhaseFunction(const Vector& v1, const Vector& v2, Float g) const {
	    Float costheta = dot(-v1, v2);
        Float gSq = g*g;
        Float num = 1.0 - gSq;
        Float den = std::pow(1.0 + gSq - 2.0 * g * costheta, 1.5);

    	return 0.5 * (num / den);
    }
    
	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
		if (m_ready)
			return true;

		if (!scene->getIntegrator()->getClass()
				->derivesFrom(MTS_CLASS(SampleIntegrator))) {
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
				index = (int) i;
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

        if (m_dumpIrrtree && m_dumpIrrtreePath.length() > 0) {
            Log(EInfo, "Starting to dump irradiance tree to %s", m_dumpIrrtreePath.c_str());
            m_octree->dumpOBJ(m_dumpIrrtreePath);
            Log(EInfo, "Dump finished");
        }

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
	Spectrum m_mfp, m_sigmaTr, m_zr, m_zv, m_alphaPrime;
	Spectrum m_sigmaSPrime, m_sigmaTPrime, m_D, m_ssFactor;
    Spectrum m_sigmaT, m_invSigmaT, m_negSigmaT;
    Float m_invSigmaTMin;
    bool m_useMartelliD;
	ref<IrradianceOctree> m_octree;
	ref<ParallelProcess> m_proc;
	int m_octreeResID, m_octreeIndex;
	int m_maxDepth;
	int m_irrSamples;
	bool m_irrIndirect;
	bool m_ready, m_requireSample, m_singleScattering, m_dumpIrrtree;
    std::string m_dumpIrrtreePath;
    mutable ThreadLocal<Random> m_random;
    bool m_useTextures;
    ref<Texture> m_zvTex;
    ref<Texture> m_zrTex;
    ref<Texture> m_sigmaTrTex;
    ref<Bitmap> m_zrBitmap;
    ref<Bitmap> m_sigmaTrBitmap;
    Float m_texUScaling;
    Float m_texVScaling;
    /* Indicates if a look-up-table should be created and used for Rd */
    bool m_useRdLookUpTable;
    /* Look-up-table for Rd, indexed by the distance r. */
    ref<LUTType> m_RdLookUpTable;
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

MTS_IMPLEMENT_CLASS_S(IsotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicDipole, "Isotropic dipole model");
MTS_NAMESPACE_END
