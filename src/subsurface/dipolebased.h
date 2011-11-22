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

#if !defined(__DIPOLE_BASED_H)
#define __DIPOLE_BASED_H

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <boost/bind.hpp>
#include <boost/timer.hpp>
#include "../medium/materials.h"
#include "irrtree.h"

//#define NO_SSE_QUERY

MTS_NAMESPACE_BEGIN

/* Relative bound on what is still accepted as roundoff 
   error -- be quite tolerant */
#if defined(SINGLE_PRECISION)
	#define ERROR_REQ 1e-2f
#else
	#define ERROR_REQ 1e-5
#endif


typedef SubsurfaceMaterialManager::LUTType LUTType;
typedef SubsurfaceMaterialManager::LUTRecord LUTRecord;

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *	 Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 */
template <typename Derived> class DipoleBasedSubsurface : public Subsurface {
protected:
	typedef boost::function<Spectrum (const Spectrum)> EvaluationFunction;

public:
	DipoleBasedSubsurface(const Properties &props) 
		: Subsurface(props) {
		irrOctreeMutex->lock();
		m_octreeIndex = irrOctreeIndex++;
		irrOctreeMutex->unlock();
		
		/* How many samples should be taken when estimating 
		   the irradiance at a given point in the scene? */
		m_irrSamples = props.getInteger("irrSamples", 32);
				
		/* When estimating the irradiance at a given point, should indirect illumination be included
		   in the final estimate? */
		m_irrIndirect = props.getBoolean("irrIndirect", true);

		/* Multiplicative factor, which can be used to adjust the number of
		   irradiance samples */
		m_sampleMultiplier = props.getFloat("sampleMultiplier", 2.0f);
		/* Error threshold - lower means better quality */
		m_minDelta = props.getFloat("quality", 0.1f);
		/* Max. depth of the created octree */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Should the irrtree be dumped? */
		m_dumpIrrtree = props.getBoolean("dumpIrrtree", false);
		m_dumpIrrtreePath = props.getString("dumpIrrtreePath", "");
		/* Multiplicative factor for the subsurface term */
		m_ssFactor = props.getSpectrum("ssFactor", Spectrum(1.0f));
		/* Asymmetry parameter of the phase function */
		m_g = props.getFloat("g", 0);
		/* alternative diffusion coefficient */
		m_useMartelliD = props.getBoolean("useMartelliDC", true);
		/* look-up-table */
		m_useLookUpTable = props.getBoolean("useLookUpTable", false);
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

		lookupMaterial(props, m_sigmaS, m_sigmaA, &m_eta);
	}

	DipoleBasedSubsurface(Stream *stream, InstanceManager *manager) 
	 : Subsurface(stream, manager) {
		m_sigmaS = Spectrum(stream);
		m_sigmaA = Spectrum(stream);
		m_ssFactor = Spectrum(stream);
		m_g = stream->readFloat();
		m_eta = stream->readFloat();
		m_sampleMultiplier = stream->readFloat();
		m_minDelta = stream->readFloat();
		m_maxDepth = stream->readInt();
		m_octreeIndex = stream->readInt();
		m_irrSamples = stream->readInt();
		m_irrIndirect = stream->readBool();
		m_useMartelliD = stream->readBool();
		m_useLookUpTable = stream->readBool();
		m_lutResolution = stream->readFloat();
		m_errThreshold = stream->readFloat();
		m_mcIterations = stream->readInt();

		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~DipoleBasedSubsurface() {
		if (m_octreeResID != -1)
			Scheduler::getInstance()->unregisterResource(m_octreeResID);
	}

	/// Cast to the derived class
	inline Derived *cast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *cast() const {
		return static_cast<const Derived *>(this);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		if (m_octreeResID != -1)
			proc->bindResource(formatString("irrOctree%i", m_octreeIndex), m_octreeResID);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		m_sigmaS.serialize(stream);
		m_sigmaA.serialize(stream);
		m_ssFactor.serialize(stream);
		stream->writeFloat(m_g);
		stream->writeFloat(m_eta);
		stream->writeFloat(m_sampleMultiplier);
		stream->writeFloat(m_minDelta);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_octreeIndex);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
		stream->writeBool(m_useMartelliD);
		stream->writeBool(m_useLookUpTable);
		stream->writeFloat(m_errThreshold);
		stream->writeFloat(m_lutResolution);
		stream->writeInt(m_mcIterations);
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		if (!m_ready || m_ssFactor.isZero())
			return Spectrum(0.0f);

		Spectrum Mo;
		if (m_useLookUpTable)
			Mo = cast()->MoLUT(scene, sampler, its, d, depth);
		else
			Mo = cast()->Mo(scene, sampler, its, d, depth);

		Spectrum Lo;
		if (m_eta == 1.0f) {
			Lo = Mo * m_ssFactor * INV_PI;
		} else {
			const Normal &n = its.shFrame.n;
			Float Ft = 1.0f - fresnel(absDot(n, d));
			Lo = Mo * m_ssFactor * INV_PI * (Ft / m_Fdr);
		}

		return Lo;
	}

	void configure() {
		boost::timer timer;
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

		/* Average transmittance at the boundary */
		m_Fdt = 1.0f - m_Fdr;

		/* Approximate dipole boundary condition term */
		m_A = (1 + m_Fdr) / m_Fdt;

		if (m_eta == 1.0f) {
			m_Fdr = (Float) 0.0f;
			m_Fdt = (Float) 1.0f;
		}

		/* Reduced albedo */
		m_alphaPrime = m_sigmaSPrime / m_sigmaTPrime;

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

		/* let the derived class configure the dipoles */
		cast()->configureDipoles();
		/* let the derived class configure a LUT */
		if (m_useLookUpTable)
			cast()->configureLUT();
	}

	/* Unpolarized fresnel reflection term for dielectric materials
	 */
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

	/* The reflectance at a particular surface position results from the
	 * irradiance from everywhere on the surface. However, only points
	 * with a maximum distance rMax will have significant impact on the
	 * result. With the help of Monte-Carlo integration this method trys
	 * to find this distance.
	 */
	Float findRMax() {
		boost::timer timer;
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
			Rd_A += cast()->getRd(r);
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
				Rd_APrime += cast()->getRd(r);
			}
			Float APrime = 4 * rMax * rMax;
			Rd_APrime = APrime * Rd_APrime * m_alphaPrime * inv4Pi / (Float)(numSamples - 1);
			err = (Rd_A - Rd_APrime) / Rd_A;
		}
		Log(EDebug, "Maximum distance for sampling surface is %f with an error of %f (took %.0fs)",
			rMax, m_errThreshold, timer.elapsed());
		return rMax;
	}

	/* Creates a look-up-table (for dMoR) that is indexed by the
	 * squared distance to the sample. For the resolution and rMax
	 * parameters the corresponding member fields are used. If rMax
	 * is was not defined by the user, it is now calculated.
	 */
	ref<LUTType> createLUT(const EvaluationFunction &f) {
		if (!m_rMaxPredefined) {
			m_rMax = findRMax();
			m_rMaxPredefined = true;
		}

		return createLUT(f, m_lutResolution, m_rMax);
	}

	/* Creates a look-up-table (for dMoR) that is indexed by the
	 * squared distance to the sample. The resolution, that is the
	 * the distance between two samples, and the maximum sampling
	 * distance is required.
	 */
	ref<LUTType> createLUT(const EvaluationFunction &f, Float resolution, Float rMax) {
		boost::timer timer;
		Log(EDebug, "Building LUT for dipole based subsurafce scattering..");

		/* build rho_dt look-up-table */
		const int numEntries = (int) (rMax / resolution) + 1;
		ref<LUTType> lut = new LUTType(numEntries);

		for (int i=0; i<numEntries; ++i) {
			const Spectrum r(i * resolution);
			lut->at(i) = f(r);
		}

		Log(EDebug, "Created a look-up-table with %i entries (took %.2fs)", numEntries, timer.elapsed());
		return lut;
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

protected:
	Float m_minMFP, m_sampleMultiplier;
	Float m_Fdr, m_Fdt, m_A, m_minDelta, m_g, m_eta;
	Spectrum m_sigmaS, m_sigmaA;
	Spectrum m_mfp, m_sigmaTr, m_alphaPrime;
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
	bool m_ready, m_dumpIrrtree;
	std::string m_dumpIrrtreePath;
	/* Indicates if a look-up-table should be created and used for Rd */
	bool m_useLookUpTable;
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

MTS_NAMESPACE_END

#endif /* __DIPOLE_BASED_H */
