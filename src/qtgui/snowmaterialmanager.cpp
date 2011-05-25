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
        // if not found, add new
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
            properties.setBoolean("addMultipleScattering", srs.hkUseMultipleScattering);
            properties.setFloat("ssFactor", srs.hkSingleScatteringFactor);
            properties.setFloat("drFactor", srs.hkMultipleScatteringFactor);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        }

        if (subsurfaceMode == ENoSubSurface) {
            subsurface = NULL; 
        } else if (subsurfaceMode == EJensenDipoleBSSRDF) {
            properties.setPluginName("dipole");
            properties.setSpectrum("ssFactor", Spectrum(srs.dipoleDensityFactor));
            properties.setFloat("sampleMultiplier", srs.dipoleSampleFactor);
            properties.setBoolean("addSingleScattering", srs.dipoleUseSingleScattering);
            properties.setBoolean("useMartelliD", srs.dipoleMartelliDC);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (subsurfaceMode == EJensenMultipoleBSSRDF) {
            properties.setPluginName("multipole");
            properties.setSpectrum("ssFactor", Spectrum(srs.multipoleDensityFactor));
            properties.setFloat("sampleMultiplier", srs.multipoleSampleFactor);
            properties.setBoolean("addSingleScattering", srs.dipoleUseSingleScattering);
            properties.setFloat("slabThickness", srs.multipoleSlabThickness);
            properties.setInteger("extraDipoles", srs.multipoleExtraDipoles);
            properties.setBoolean("useMartelliD", srs.multipoleMartelliDC);
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

bool SnowMaterialManager::isMadeOfSnow(Shape * shape) {
    ShapeMap::iterator it = snowShapes.find(shape);
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
            Shape *s = it->first;
            ShapeEntry &entry = it->second;
            if (entry.madeOfSnow && s != NULL) {
                BSDF* bsdf = s->getBSDF();
                Subsurface* subsurface = s->getSubsurface();

                oss << "  " << s->getName() << ":" << std::endl
                << "    BSDF: " << std::endl << (bsdf == NULL ? "None" : bsdf->toString()) << std::endl
                << "    Subsurface: " << std::endl << (subsurface == NULL ? "None" : subsurface->toString()) << std::endl;
            }
        }
		oss	<< "]";
		return oss.str();
}

MTS_NAMESPACE_END
