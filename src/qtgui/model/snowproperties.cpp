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
        Float _ior, Float _g)
     :  grainsize(_grainsize), density(_density), ior(_ior), g(_g) {
    configure();
}

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
    lastPreset = EFreshNewSnow;
    grainsize = 0.05f / 1000.0f; // mm -> m
    density = 70;
    ior = 1.32;
    g = 0.5;
    configure();
}

void SnowProperties::loadDryOlderSnowPreset() {
    lastPreset = EDryOlderSnow;
    grainsize = 0.25f / 1000.0f; // mm -> m
    density = 300;
    ior = 1.32;
    g = 0.5;
    configure();
}

void SnowProperties::loadWetOldSnowPreset() {
    lastPreset = EWetOldSnow;
    grainsize = 1.0f / 1000.0f; // mm -> m
    density = 450;
    ior = 1.32;
    g = 0.5;
    configure();
}

void SnowProperties::configure() {
    sigmaA = getSigmaA(iceSigmaA, density, iceDensity);
    sigmaT = getAsymptoticExtCoeff(sigmaA, grainsize, density, iceDensity);
    sigmaS = sigmaT - sigmaA;
    singleScatteringAlbedo = sigmaS / sigmaT;
}

std::string SnowProperties::toString() {
		std::ostringstream oss;
		oss << "SnowProperties[" << std::endl
            << "  sigmaA = " << sigmaA.toString() << std::endl
            << "  sigmaS = " << sigmaS.toString() << std::endl
            << "  sigmaT = " << sigmaT.toString() << std::endl
            << "  sigmaA (ice) = " << iceSigmaA.toString() << std::endl
            << "  ss. Albedo = " << singleScatteringAlbedo.toString() << std::endl
		    << "]";
		return oss.str();
}

MTS_IMPLEMENT_CLASS(SnowProperties, false, Object)
MTS_NAMESPACE_END
