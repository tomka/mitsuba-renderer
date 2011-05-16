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
        // try to find shape in backup lists
        BSDFMap::iterator bsdfEntry = originalBSDFs.find(shape);
        SubsurfaceMap::iterator subsurfaceEntry = originalSubsurfaces.find(shape);

        // if not found, add them
        if (bsdfEntry == originalBSDFs.end()) {
            originalBSDFs[shape] = bsdfOld;
            if (bsdfOld != NULL)
                bsdfOld->incRef();
        }
        if (subsurfaceEntry == originalSubsurfaces.end()) {
            originalSubsurfaces[shape] = subsurfaceOld;
            if (subsurfaceOld != NULL)
                subsurfaceOld->incRef();
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
        if (surfaceMode == ENoSurface) {
            bsdf = NULL;
        } else if (surfaceMode == EWiscombeWarrenAlbedo) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth", 2.0f); // ToDo: Make dynamic
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (surfaceMode == EWiscombeWarrenBRDF) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth", 2.0f); // ToDo: Make dynamic
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (surfaceMode == EHanrahanKruegerBRDF) {
            properties.setPluginName("hanrahankrueger");
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        }

        if (subsurfaceMode == ENoSubSurface) {
            subsurface = NULL; 
        } else if (subsurfaceMode == EJensenDipoleBSSRDF) {
            properties.setPluginName("dipole");
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (subsurfaceMode == EJensenMultipoleBSSRDF) {
            properties.setPluginName("multipole");
            properties.setFloat("slabThickness", 0.1); // ToDo: Make dynamic
            properties.setInteger("extraDipoles", context->multipoleDipoles);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (subsurfaceMode == EJakobADipoleBSSRDF) {
            properties.setPluginName("adipole");
            properties.setString("D", getFlakeDistribution());
            properties.setFloat("sigmaTn", 1.0f);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        }

        /* initialize new materials */
        if (bsdf) {
            bsdf->setParent(shape);
            bsdf->configure();
        }
        if (subsurface) {
            subsurface->setParent(shape);
            subsurface->configure();
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
        materialMap[shape] = std::pair<BSDF*, Subsurface*>(bsdf, subsurface);
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
        // try to find shape in backup lists
        BSDFMap::iterator bsdfEntry = originalBSDFs.find(shape);
        SubsurfaceMap::iterator subsurfaceEntry = originalSubsurfaces.find(shape);
        MaterialMap::iterator materialEntry = materialMap.find(shape);

        // if found, use materials
        if (bsdfEntry != originalBSDFs.end()) {
            shape->setBSDF( (*bsdfEntry).second );
        }
        if (subsurfaceEntry != originalSubsurfaces.end()) {
            Subsurface *old_ss = shape->getSubsurface();
            Subsurface *subsurface = (*subsurfaceEntry).second;
            shape->setSubsurface(subsurface);
            if (subsurface != NULL)
                context->scene->addSubsurface(subsurface);
            if (old_ss != NULL)
                context->scene->removeSubsurface(old_ss);
        }
        if (materialEntry != materialMap.end()) {
            materialMap.erase(materialEntry);
        }
        setMadeOfSnow(shape, false);
        std::cerr << "[Snow Material Manager] Reset material on shape " << shape->getName() << std::endl;
}

bool SnowMaterialManager::isMadeOfSnow(Shape * shape) {
    if (snowShapes.find(shape) != snowShapes.end())
        return snowShapes[shape];
    else
        return false;
}

void SnowMaterialManager::setMadeOfSnow(Shape *shape, bool snow) {
    snowShapes[shape] = snow;
}

std::string SnowMaterialManager::toString() {
		std::ostringstream oss;
		oss << "SnowMaterialManager[" << std::endl;

        for (ShapeMap::iterator it = snowShapes.begin(); it != snowShapes.end(); it++) {
            Shape *s = it->first;
            if (it->second) {
                BSDF* bsdf = materialMap[s].first;
                Subsurface* subsurface = materialMap[s].second;

                oss << "  " << s->getName() << ":" << std::endl
                << "    BSDF: " << std::endl << (bsdf == NULL ? "None" : bsdf->toString()) << std::endl
                << "    Subsurface: " << std::endl << (subsurface == NULL ? "None" : subsurface->toString()) << std::endl;
            }
        }

		oss	<< "]";
		return oss.str();
}

MTS_NAMESPACE_END
