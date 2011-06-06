#ifndef SNOWMATERALMANAGER_H
#define SNOWMATERALMANAGER_H

#include "common.h"
#include <mitsuba/mitsuba.h>
#include "model/snowproperties.h"
#include <map>

MTS_NAMESPACE_BEGIN

class SnowMaterialManager {

    struct ShapeEntry {
        bool madeOfSnow;
        BSDF *originalBSDF;
        Subsurface *originalSubsurface;

        ShapeEntry() : madeOfSnow(false),
            originalBSDF(NULL), originalSubsurface(NULL) { }
    };

    typedef std::map<Shape*, ShapeEntry> ShapeMap;

    // a toggle for every shape if it has a snow material
    ShapeMap snowShapes;

    /* A storage point for a snow diffusion profile. This is mainly used for
     * realtime SSS. */
    ref<Bitmap> diffusionProfileCache;
    Float diffusionProfileRmax;

public:
    /**
     * Add a material defined by 'mode' to the shape. The required snow
     * properties are passed, too.
     */
    void replaceMaterial(Shape *shape, SceneContext *context);

    /**
     * Reset material (BSDF or Subsurface) of the shape to its orignal state.
     */
    void resetMaterial(Shape *shape, SceneContext *context);

    /**
     * Query if a particular shape is made of snow.
     */
    bool isMadeOfSnow(Shape *shape);

    /**
     * Removes all links to a shape.
     */
    void removeShape(Shape *shape);

    /**
     * Get a string description of the snow material manager.
     */
    std::string toString();

    /**
     * Get last diffusion profile calculated.
     */
    std::pair< ref<Bitmap>, Float > getCachedDiffusionProfile() const;

    /**
     * Check if there is a diffusion profile available that is ready for usage.
     */
    bool hasCachedDiffusionProfile() const;

    /**
     * Replace the currently cached diffusion profile with a new one.
     */
    void refreshDiffusionProfile(const SceneContext *context);

    /**
     * Calculate the dipole diffusion napproximation.
     */
    Spectrum getRd(Spectrum &r, Spectrum &sigmaTr, Spectrum &zv, Spectrum &zr);

protected:
    /**
     * Enable/Disable snow material for a specific shape.
     */
    void setMadeOfSnow(Shape *shape, bool snow);

    /**
     * Get different flake distributions for the anisotropic BSSRDF.
     */
    std::string getFlakeDistribution();
};

MTS_NAMESPACE_END

#endif /* SNOWMATERALMANAGER_H */
