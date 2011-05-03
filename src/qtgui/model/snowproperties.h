#ifndef SNOWPROPERTIES_H
#define SNOWPROPERTIES_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

struct SnowProperties {
    MTS_DECLARE_CLASS()

    /* Some preset configurations */
    enum EPreset {
        EFreshNewSnow = 0,
        EDryOlderSnow,
        EWetOldSnow  
    };

    /* grain diameter in m */
    Float grainsize;
    /* density in kg / m^3 */
    Float density;
    /* IOR */
    Float ior;
    /* asymmetry factor g (mean cosine of phase function */
    Float g;
    /* absorption coefficient */
    Spectrum sigmaA;
    /* scattering coefficient */
    Spectrum sigmaS;
    /* extinction coefficient, sum of sigmaA and sigmaS */
    Spectrum sigmaT;
    /* absorbtion coefficient of ice */
    static Spectrum iceSigmaA;
    /* density of ice */
    static Float iceDensity;

    SnowProperties();

    SnowProperties(EPreset preset);

    SnowProperties(Float _grainsize, Float _density, Float _ior,
            Float _g, Float _sigmaA, Float _sigmaS);

    void loadPreset(EPreset preset);

    void loadFreshNewSnowPreset();

    void loadDryOlderSnowPreset();

    void loadWetOldSnowPreset();

    void loadCoefficients();
};

MTS_NAMESPACE_END

#endif // SNOWPROPERTIES_H
