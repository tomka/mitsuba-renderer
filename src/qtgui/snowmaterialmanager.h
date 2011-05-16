#ifndef SNOWMATERALMANAGER_H
#define SNOWMATERALMANAGER_H

#include "common.h"
#include <mitsuba/mitsuba.h>
#include "model/snowproperties.h"
#include <map>

MTS_NAMESPACE_BEGIN

class SnowMaterialManager {
    typedef std::map<Shape*, BSDF*> BSDFMap;
    typedef std::map<Shape*, Subsurface*> SubsurfaceMap;
    typedef std::map<Shape*, bool> ShapeMap;
    typedef std::map<Shape*, std::pair<BSDF*, Subsurface*> > MaterialMap; 

    // original BSDFs of shapes
    BSDFMap originalBSDFs;
    // original Subsurface integrators of shapes
    SubsurfaceMap originalSubsurfaces;
    // a toggle for every shape if it has a snow material
    ShapeMap snowShapes;
    // a map of the currently assigned materials
    MaterialMap materialMap;

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
    bool isMadeOfSnow(Shape *);

    /**
     * Get a string description of the snow material manager.
     */
    std::string toString();

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
