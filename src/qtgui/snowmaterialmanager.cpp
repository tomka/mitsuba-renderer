#include "snowmaterialmanager.h"
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

void SnowMaterialManager::replaceMaterial(Shape *shape, SceneContext *context) {
        BSDF *bsdf = shape->getBSDF(); 
        Subsurface *subsurface = shape->getSubsurface(); 
        // try to find shape in backup lists
        BSDFMap::iterator bsdfEntry = originalBSDFs.find(shape);
        SubsurfaceMap::iterator subsurfaceEntry = originalSubsurfaces.find(shape);

        // if not found, add them
        if (bsdfEntry == originalBSDFs.end()) {
            originalBSDFs[shape] = bsdf;
            if (bsdf != NULL)
                bsdf->incRef();
        }
        if (subsurfaceEntry == originalSubsurfaces.end()) {
            originalSubsurfaces[shape] = subsurface;
            if (subsurface != NULL)
                subsurface->incRef();
        }

        // remove all children that are either subsurface materials or BSDFs
        bsdf = NULL; 
        subsurface = NULL; 
        
        PluginManager *pluginManager = PluginManager::getInstance();
        ERenderMode mode = context->snowRenderMode;
        SnowProperties &snow = context->snow;
        Properties properties;
        properties.setFloat("g", snow.g);
        properties.setSpectrum("sigmaA", snow.sigmaA);
        properties.setSpectrum("sigmaS", snow.sigmaS);
        properties.setSpectrum("sigmaT", snow.sigmaT);
        properties.setFloat("eta", snow.ior); // ToDo: eta is actually the relative IOR (no prob w/ air)
        if (mode == EWiscombeWarrenAlbedo) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth", 2.0f); // ToDo: Make dynamic
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (mode == EWiscombeWarrenBRDF) {
            properties.setPluginName("wiscombe");
            properties.setFloat("depth", 2.0f); // ToDo: Make dynamic
            properties.setSpectrum("singleScatteringAlbedo", snow.singleScatteringAlbedo);
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (mode == EHanrahanKruegerBRDF) {
            properties.setPluginName("hanrahankrueger");
            bsdf = static_cast<BSDF *> (pluginManager->createObject(
                BSDF::m_theClass, properties));
        } else if (mode == EJensenBSSRDF) {
            properties.setPluginName("dipole");
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        } else if (mode == EJensenMultipoleBSSRDF) {
            properties.setPluginName("multipole");
            properties.setFloat("slabThickness", 0.1); // ToDo: Make dynamic
            properties.setInteger("extraDipoles", context->multipoleDipoles);
            subsurface = static_cast<Subsurface *> (pluginManager->createObject(
                Subsurface::m_theClass, properties));
        }

        if (bsdf) {
            bsdf->setParent(shape);
            bsdf->configure();
        }
        if (subsurface) {
            subsurface->setParent(shape);
            subsurface->configure();
        }

        shape->setBSDF(bsdf);
        shape->setSubsurface(subsurface);
        // allow the shape to react to this changes
        shape->configure();

        // if a subsurface material has been selected, inform the scene about it
        if (subsurface) {
            context->scene->addSubsurface(subsurface);
            /* if the subsurface integrator previously used (if any) is not
             * needed by other shapes, we can remove it for now.
             */
            Subsurface *old_ss = originalSubsurfaces[shape];
            if (old_ss != NULL ) {
                context->scene->removeSubsurface(old_ss);
            }
        }
        setMadeOfSnow(shape, true);
        std::string bsdfName = (bsdf == NULL) ? "None" : bsdf->getClass()->getName();
        std::string subsurfaceName = (subsurface == NULL) ? "None" : subsurface->getClass()->getName();
        std::cerr << "[Snow Material Manager] Replaced material of shape \"" << shape->getName() << "\"" << std::endl
                  << "\tnew BSDF: " << bsdfName << std::endl
                  << "\tnew Subsurface: " << subsurfaceName << std::endl;
}

void SnowMaterialManager::resetMaterial(Shape *shape, SceneContext *context) {
        // try to find shape in backup lists
        BSDFMap::iterator bsdfEntry = originalBSDFs.find(shape);
        SubsurfaceMap::iterator subsurfaceEntry = originalSubsurfaces.find(shape);

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

MTS_NAMESPACE_END
