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
    // original BSDFs of shapes
    BSDFMap originalBSDFs;
    // original Subsurface integrators of shapes
    SubsurfaceMap originalSubsurfaces;
    // a toggle for every shape if it has a snow material
    std::map<Shape *, bool> snowShapes;

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

protected:
    /**
     * Enable/Disable snow material for a specific shape.
     */
    void setMadeOfSnow(Shape *shape, bool snow);
};

MTS_NAMESPACE_END

#endif /* SNOWMATERALMANAGER_H */
