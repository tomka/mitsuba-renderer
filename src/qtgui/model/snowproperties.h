#ifndef SNOWPROPERTIES_H
#define SNOWPROPERTIES_H

#include <mitsuba/mitsuba.h>
#include <math/snowmath.h>

MTS_NAMESPACE_BEGIN

struct SnowProperties {
    MTS_DECLARE_CLASS()

    /* Some preset configurations */
    enum EPreset {
        EFreshNewSnow = 0,
        EDryOlderSnow,
        EWetOldSnow  
    };

    /* grain diameter in um */
    Float grainsize;
    /* density in kg / m^3 */
    Float density;
    /* IOR */
    Float ior;
    /* asymetry factor g (mean cosine of phase function */
    Float g;
    /* absorption coefficient */
    Spectrum sigmaA;
    /* scattering coefficient */
    Spectrum sigmaS;

    SnowProperties();

    SnowProperties(EPreset preset);

    SnowProperties(Float _grainsize, Float _density, Float _ior,
            Float _g, Float _sigmaA, Float _sigmaS);

    void loadPreset(EPreset preset);

    void loadFreshNewSnowPreset();

    void loadDryOlderSnowPreset();

    void loadWetOldSnowPreset();
};

MTS_NAMESPACE_END

#endif // SNOWPROPERTIES_H
