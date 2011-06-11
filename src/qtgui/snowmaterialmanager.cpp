#include "common.h"
#include "snowmaterialmanager.h"
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

void SnowMaterialManager::replaceMaterial(Shape *shape, SceneContext *context) {
        // get currently selected material
        BSDF *bsdfOld = shape->getBSDF();
        Subsurface *subsurfaceOld = shape->getSubsurface();

        // try to find shape in store 
        ShapeMap::iterator it = snowShapes.find(shape);
        /* If not found, add new.
         */
        if (it == snowShapes.end()) {
            ShapeEntry newEntry;

            newEntry.originalBSDF = bsdfOld;
            if (bsdfOld != NULL)
                bsdfOld->incRef();

            newEntry.originalSubsurface = subsurfaceOld;
            if (subsurfaceOld != NULL)
                subsurfaceOld->incRef();

            shape->incRef();
            snowShapes[shape] = newEntry;
        }

        PluginManager *pluginManager = PluginManager::getInstance();
        SnowProperties &snow = context->snow;
        ESurfaceRenderMode surfaceMode = context->snowRenderSettings.surfaceRenderMode;
        ESubSurfaceRenderMode subsurfaceMode = context->snowRenderSettings.subsurfaceRenderMode;

        // common properties
        Properties properties;
        properties.setFloat("g", snow.g);
        properties.setSpectrum("sigmaA", snow.sigmaA);
        properties.setSpectrum("sigmaS", snow.sigmaS);
        properties.setSpectrum("sigmaT", snow.sigmaT);
        properties.setFloat("eta", snow.ior); // ToDo: eta is actually the relative IOR (no prob w/ air)
 
        BSDF *bsdf = NULL;
        Subsurface *subsurface = NULL;
        SnowRenderSettings &srs = context->snowRenderSettings;

        if (surfaceMode == ENoSurface) {
            bsdf = NULL;
        } else if (surfaceMode == EWiscombeWarrenAlbedo) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth",srs.wiscombeDepth);
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (surfaceMode == EWiscombeWarrenBRDF) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth", srs.wiscombeDepth);
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (surfaceMode == EHanrahanKruegerBRDF) {
            properties.setPluginName("hanrahankrueger");
            properties.setFloat("ssFactor", srs.hkSingleScatteringFactor);
            properties.setFloat("drFactor", srs.hkMultipleScatteringFactor);
            properties.setBoolean("diffuseReflectance", srs.hkUseMultipleScattering);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        }

        if (subsurfaceMode == ENoSubSurface) {
            subsurface = NULL; 
        } else if (subsurfaceMode == EJensenDipoleBSSRDF) {
            properties.setPluginName("dipole");
            properties.setSpectrum("ssFactor", Spectrum(srs.dipoleDensityFactor));
            properties.setFloat("sampleMultiplier", srs.dipoleSampleFactor);
            properties.setBoolean("singleScattering", srs.dipoleUseSingleScattering);
            properties.setBoolean("useMartelliD", srs.dipoleMartelliDC);
            properties.setBoolean("useTexture", srs.dipoleTexture);
            properties.setBoolean("dumpIrrtree", srs.dipoleDumpIrrtree);
            properties.setString("dumpIrrtreePath", srs.dipoleDumpIrrtreePath);
            if (srs.dipoleLutPredefineRmax) {
                properties.setFloat("lutRmax", srs.dipoleLutRmax);
            } else {
                properties.setInteger("mcIterations", srs.dipoleLutMCIterations);
            }
            if (srs.dipoleTexture) {
                properties.setString("zrFilename", srs.dipoleZrTexture);
                properties.setString("sigmaTrFilename", srs.dipoleSigmaTrTexture);
                properties.setFloat("texUScaling", srs.dipoleTextureUScaling);
                properties.setFloat("texVScaling", srs.dipoleTextureVScaling);
            }
            properties.setBoolean("useLookUpTable", srs.dipoleUseLut);
            properties.setFloat("lutResolution", srs.dipoleLutResolution);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (subsurfaceMode == EJensenMultipoleBSSRDF) {
            properties.setPluginName("multipole");
            properties.setSpectrum("ssFactor", Spectrum(srs.multipoleDensityFactor));
            properties.setFloat("sampleMultiplier", srs.multipoleSampleFactor);
            properties.setBoolean("singleScattering", srs.dipoleUseSingleScattering);
            properties.setFloat("slabThickness", srs.multipoleSlabThickness);
            properties.setInteger("extraDipoles", srs.multipoleExtraDipoles);
            properties.setBoolean("useMartelliD", srs.multipoleMartelliDC);
            properties.setBoolean("useLookUpTable", srs.multipoleUseLut);
            properties.setFloat("lutResolution", srs.multipoleLutResolution);
            if (srs.multipoleLutPredefineRmax) {
                properties.setFloat("lutRmax", srs.multipoleLutRmax);
            } else {
                properties.setInteger("mcIterations", srs.multipoleLutMCIterations);
            }
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (subsurfaceMode == EJakobADipoleBSSRDF) {
            properties.setPluginName("adipole");
            properties.setSpectrum("ssFactor", Spectrum(srs.adipoleDensityFactor));
            properties.setFloat("sampleMultiplier", srs.adipoleSampleFactor);
            QString D = QString::fromStdString(srs.adipoleD);
            if (D.trimmed().length() == 0)
                properties.setString("D", getFlakeDistribution());
            else
                properties.setString("D", srs.adipoleD);
            properties.setFloat("sigmaTn", srs.adipoleSigmaTn);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        }

        /* initialize new materials */
        if (bsdf) {
            bsdf->setParent(shape);
            bsdf->configure();
            bsdf->incRef();
        }
        if (subsurface) {
            subsurface->setParent(shape);
            subsurface->configure();
            subsurface->incRef();
            // if a subsurface material has been selected, inform the scene about it
            context->scene->addSubsurface(subsurface);
        }

        shape->setBSDF(bsdf);
        shape->setSubsurface(subsurface);
        // allow the shape to react to this changes
        shape->configure();

        /* if the subsurface integrator previously used (if any) is not
         * needed by other shapes, we can remove it for now.
         */
        if (subsurfaceOld != NULL ) {
            context->scene->removeSubsurface(subsurfaceOld);
        }

        setMadeOfSnow(shape, true);
        std::string bsdfName = (bsdf == NULL) ? "None" : bsdf->getClass()->getName();
        std::string subsurfaceName = (subsurface == NULL) ? "None" : subsurface->getClass()->getName();
        std::cerr << "[Snow Material Manager] Replaced material of shape \"" << shape->getName() << "\"" << std::endl
                  << "\tnew BSDF: " << bsdfName << std::endl
                  << "\tnew Subsurface: " << subsurfaceName << std::endl;
}

std::string SnowMaterialManager::getFlakeDistribution() {
    /* clamped cosin^20 flake distribution */
    // return "0.01314, -0.00014, 0.00061, -0.00014, 0.01295, -0.00018, 0.00061, -0.00018, -0.07397";
    /* sine^20 flake distribution */
    return "1.6307, -0.00049, 0.00069, -0.00049, 1.63148, 0.00001, 0.00067, 0.00002, 2.12596";
}

void SnowMaterialManager::resetMaterial(Shape *shape, SceneContext *context) {
        ShapeMap::iterator it = snowShapes.find(shape);
        if (it == snowShapes.end()) {
            SLog(EWarn, "Did not find requested shape to reset material.");
            return;
        }

        setMadeOfSnow(shape, false);
        ShapeEntry &e = it->second;

        // if found, use materials
        BSDF *bsdf = e.originalBSDF;
        if (bsdf != NULL) {
            bsdf->setParent(shape);
            bsdf->configure();
        }
        shape->setBSDF( bsdf );

        Subsurface *old_ss = shape->getSubsurface();
        Subsurface *subsurface = e.originalSubsurface;

        if (subsurface != NULL) {
            subsurface->setParent(shape);
            subsurface->configure();
            context->scene->addSubsurface(subsurface);
        }
        if (old_ss != NULL)
            context->scene->removeSubsurface(old_ss);

        shape->setSubsurface(subsurface);
        // allow the shape to react to this changes
        shape->configure();

        std::cerr << "[Snow Material Manager] Reset material on shape " << shape->getName() << std::endl;
}

bool SnowMaterialManager::isMadeOfSnow(const Shape * shape) const {
    ShapeMap::const_iterator it = snowShapes.find(shape);
    if (it != snowShapes.end())
        return it->second.madeOfSnow;
    else
        return false;
}

void SnowMaterialManager::removeShape(Shape *shape) {
    ShapeMap::iterator it = snowShapes.find(shape);
    if (it != snowShapes.end())
        snowShapes.erase(it);
}

void SnowMaterialManager::setMadeOfSnow(Shape *shape, bool snow) {
    ShapeMap::iterator it = snowShapes.find(shape);
    if (it == snowShapes.end())
        return;

    it->second.madeOfSnow = snow;
}

std::string SnowMaterialManager::toString() {
		std::ostringstream oss;
		oss << "SnowMaterialManager[" << std::endl;

        for (ShapeMap::iterator it = snowShapes.begin(); it != snowShapes.end(); it++) {
            const Shape *s = it->first;
            ShapeEntry &entry = it->second;
            if (entry.madeOfSnow && s != NULL) {
                const BSDF* bsdf = s->getBSDF();
                const Subsurface* subsurface = s->getSubsurface();

                oss << "  " << s->getName() << ":" << std::endl
                << "    BSDF: " << std::endl << (bsdf == NULL ? "None" : indent(bsdf->toString(), 3)) << std::endl
                << "    Subsurface: " << std::endl << (subsurface == NULL ? "None" : indent(subsurface->toString(), 3)) << std::endl;
            }
        }
		oss	<< "]";
		return oss.str();
}

std::pair< ref<Bitmap>, Float > SnowMaterialManager::getCachedDiffusionProfile() const {
    return std::make_pair(diffusionProfileCache, diffusionProfileRmax);
}

bool SnowMaterialManager::hasCachedDiffusionProfile() const {
    return diffusionProfileCache.get() != NULL;
}

void SnowMaterialManager::refreshDiffusionProfile(const SceneContext *context) {
    typedef SubsurfaceMaterialManager::LUTType LUTType;

    const Float errThreshold = 0.01f;
    const Float lutResolution = 0.01f;
    const SnowProperties &sp = context->snow;

    Spectrum sigmaSPrime = sp.sigmaS * (1 - sp.g);
    Spectrum sigmaTPrime = sigmaSPrime + sp.sigmaA;

    /* extinction coefficient */
    Spectrum sigmaT = sp.sigmaA + sp.sigmaS;

    /* Effective transport extinction coefficient */
    Spectrum sigmaTr = (sp.sigmaA * sigmaTPrime * 3.0f).sqrt();

    /* Reduced albedo */
    Spectrum alphaPrime = sigmaSPrime / sigmaTPrime;

    /* Mean-free path (avg. distance traveled through the medium) */
    Spectrum mfp = Spectrum(1.0f) / sigmaTPrime;

    Float Fdr;
    if (sp.ior > 1) {
        /* Average reflectance due to mismatched indices of refraction
           at the boundary - [Groenhuis et al. 1983]*/
        Fdr = -1.440f / (sp.ior * sp.ior) + 0.710f / sp.ior
            + 0.668f + 0.0636f * sp.ior;
    } else {
        /* Average reflectance due to mismatched indices of refraction
           at the boundary - [Egan et al. 1973]*/
        Fdr = -0.4399f + 0.7099f / sp.ior - 0.3319f / (sp.ior * sp.ior)
            + 0.0636f / (sp.ior * sp.ior * sp.ior);
    }

    /* Average transmittance at the boundary */
    Float Fdt = 1.0f - Fdr;

    if (sp.ior == 1.0f) {
        Fdr = (Float) 0.0f;
        Fdt = (Float) 1.0f;
    }

    /* Approximate dipole boundary condition term */
    Float A = (1 + Fdr) / Fdt;

    /* Distance of the dipole point sources to the surface */
    Spectrum zr = mfp;
    Spectrum zv = mfp * (1.0f + 4.0f/3.0f * A);

    const Spectrum invSigmaTr = 1.0f / sigmaTr;
    const Float inv4Pi = 1.0f / (4 * M_PI);
    ref<Random> random = new Random();

    /* Find Rd for the whole area by monte carlo integration. The
     * sampling area is calculated from the max. mean free path.
     * A square area around with edge length 2 * maxMFP is used
     * for this. Hene, the sampling area is 4 * maxMFP * maxMFP. */
    const int numSamples = 10000;
    Spectrum Rd_A = Spectrum(0.0f);
    for (int n = 0; n < numSamples; ++n) {
        /* do importance sampling by choosing samples distributed
         * with sigmaTr^2 * e^(-sigmaTr * r). */
        Spectrum r = invSigmaTr * -std::log( random->nextFloat() );
        Rd_A += getRd(r, sigmaTr, zv, zr);
    }
    Float Area = 4 * invSigmaTr.max() * invSigmaTr.max();
    Rd_A = Area * Rd_A * alphaPrime * inv4Pi / (Float)(numSamples - 1);
    SLog(EDebug, "After %i MC integration iterations, Rd seems to be %s", numSamples, Rd_A.toString().c_str());

    /* Since we now have Rd integrated over the whole surface we can find a valid rmax
     * for the given threshold. */
    const Float step = lutResolution;
    Float rMax = 0.0f;
    Spectrum err(std::numeric_limits<Float>::max());
    while (err.max() > errThreshold) {
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
            Rd_APrime += getRd(r, sigmaTr, zv, zr);
        }
        Float APrime = 4 * rMax * rMax;
        Rd_APrime = APrime * Rd_APrime * alphaPrime * inv4Pi / (Float)(numSamples - 1);
        err = (Rd_A - Rd_APrime) / Rd_A;
    }
    SLog(EDebug, "Maximum distance for sampling surface is %f with an error of %f", rMax, errThreshold);

    /* Create the actual look-up-table */
    const int numEntries = (int) (rMax / step) + 1;
    std::vector<Spectrum> diffusionProfile;
    for (int i=0; i<numEntries; ++i) {
        Spectrum r = Spectrum(i * step);
        diffusionProfile.push_back( getRd( r, sigmaTr, zv, zr) );
std::cerr << diffusionProfile[i].toString() << std::endl;
    }
    SLog(EDebug, "Created Rd diffusion profile with %i entries.", numEntries);
}


/// Calculate Rd based on all dipoles and the requested distance
Spectrum SnowMaterialManager::getRd(Spectrum &r, Spectrum &sigmaTr, Spectrum &zv, Spectrum &zr) { 
    const Spectrum one(1.0f);
    const Spectrum negSigmaTr = sigmaTr * (-1.0f);
    const Spectrum rSqr = r * r; 
                    
    // calulate diffuse reflectance and transmittance 
    Spectrum dr = (rSqr + zr*zr).sqrt();
    Spectrum dv = (rSqr + zv*zv).sqrt();
                    
    // the change in Rd
    Spectrum Rd =   (zr * (one + sigmaTr * dr) * (negSigmaTr * dr).exp() / (dr * dr * dr))
                  + (zv * (one + sigmaTr * dv) * (negSigmaTr * dv).exp() / (dv * dv * dv));
    return Rd;  
}

MTS_NAMESPACE_END
