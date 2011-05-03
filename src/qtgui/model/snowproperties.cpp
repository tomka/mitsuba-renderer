#include "snowproperties.h"
#include <math/snowmath.h>

MTS_NAMESPACE_BEGIN

/* init static absorbtion coefficist of ice (in 1/m) */
Spectrum SnowProperties::iceSigmaA = getSigmaAofIce();
/* init density of ice (in kg / m^-3) */
Float SnowProperties::iceDensity = 917.0f;

SnowProperties::SnowProperties() {
    loadPreset(EFreshNewSnow);
}

SnowProperties::SnowProperties(EPreset preset) {
    loadPreset(preset);
}

SnowProperties::SnowProperties(Float _grainsize, Float _density,
        Float _ior, Float _g, Float _sigmaA, Float _sigmaS)
     :  g(_g), sigmaA(_sigmaA), sigmaS(_sigmaS)
{ }

void SnowProperties::loadPreset(EPreset preset) {
    if (preset == EFreshNewSnow)
        loadFreshNewSnowPreset();
    else if (preset == EDryOlderSnow)
        loadDryOlderSnowPreset();
    else if (preset == EWetOldSnow)
        loadWetOldSnowPreset();
    else {
        Log(EWarn, "An unknown preset was requested, I'll use \"fresh new snow\" instead.");
        loadFreshNewSnowPreset();
    }
}

void SnowProperties::loadFreshNewSnowPreset() {
    grainsize = 0.05f / 1000.0f; // mm -> m
    density = 70;
    ior = 1.32;
    g = 0.874;
    loadCoefficients();
}

void SnowProperties::loadDryOlderSnowPreset() {
    grainsize = 0.25f / 1000.0f; // mm -> m
    density = 300;
    ior = 1.32;
    g = 0.874;
    loadCoefficients();
}

void SnowProperties::loadWetOldSnowPreset() {
    grainsize = 1.0f / 1000.0f; // mm -> m
    density = 450;
    ior = 1.32;
    g = 0.874;
    loadCoefficients();
}

void SnowProperties::loadCoefficients() {
    sigmaA = getSigmaA(iceSigmaA, density, iceDensity);
    sigmaT = getAsymptoticExtCoeff(sigmaA, grainsize, density, iceDensity);
    sigmaS = sigmaT - sigmaA;
}

MTS_IMPLEMENT_CLASS(SnowProperties, false, Object)
MTS_NAMESPACE_END
