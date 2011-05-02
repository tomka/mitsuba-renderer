#include "snowproperties.h"
#include <math/snowmath.h>

MTS_NAMESPACE_BEGIN

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
    grainsize = 50;
    density = 70;
    ior = 1.32;
    g = 0.874;
}

void SnowProperties::loadDryOlderSnowPreset() {
    grainsize = 250;
    density = 300;
    ior = 1.32;
    g = 0.874;
}

void SnowProperties::loadWetOldSnowPreset() {
    grainsize = 1000;
    density = 450;
    ior = 1.32;
    g = 0.874;
}

MTS_IMPLEMENT_CLASS(SnowProperties, false, Object)
MTS_NAMESPACE_END
