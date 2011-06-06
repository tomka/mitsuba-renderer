#include "snowproperties.h"
#include <math/snowmath.h>

MTS_NAMESPACE_BEGIN

/* init static absorbtion coefficist of ice (in 1/m) */
Spectrum SnowProperties::iceSigmaA = getSigmaAofIce();
/* init density of ice (in kg / m^-3) */
Float SnowProperties::iceDensity = 917.0f;

SnowProperties::SnowProperties()
    : calcMode(ELargeParticle) {
    loadPreset(EFreshNewSnow);
}

SnowProperties::SnowProperties(EPreset preset)
    : calcMode(EAsymptotic) {
    loadPreset(preset);
}

SnowProperties::SnowProperties(Float _grainsize, Float _density,
        Float _ior, Float _g)
     :  grainsize(_grainsize), density(_density), ior(_ior), g(_g), calcMode(EAsymptotic) {
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
    ior = 1.31;
    g = 0.78;
    configure();
}

void SnowProperties::loadDryOlderSnowPreset() {
    lastPreset = EDryOlderSnow;
    grainsize = 0.25f / 1000.0f; // mm -> m
    density = 300;
    ior = 1.31;
    g = 0.78;
    configure();
}

void SnowProperties::loadWetOldSnowPreset() {
    lastPreset = EWetOldSnow;
    grainsize = 1.0f / 1000.0f; // mm -> m
    density = 450;
    ior = 1.31;
    g = 0.78;
    configure();
}

void SnowProperties::configure() {
    sigmaA = getSigmaA(iceSigmaA, density, iceDensity);
    if (calcMode == EPhenomenological) {
        singleScatteringAlbedo = Spectrum(0.99f);
        //singleScatteringAlbedo = getSingleScatteringAlbedo(iceSigmaA, grainsize);
        Spectrum v0 = geRteEigenvector(singleScatteringAlbedo, g);
        sigmaT = getBarkstromExtCoeff(iceSigmaA, grainsize, density, iceDensity, v0);
        sigmaS = sigmaT - sigmaA;
    } else {
        if (calcMode == EAsymptotic) {
            sigmaT = getAsymptoticExtCoeff(sigmaA, grainsize, density, iceDensity);
        } else if (calcMode == ESnowPack) {
            Float c1 = 10.0f; // kg * m^-2
            Float c2 = 30.0f; // m^-1
            sigmaT = getSnowPackExtCoeff(density, c1, c2);
        } else if (calcMode == ELargeParticle) {
            sigmaT = getLargeParticleExtCoeff(grainsize, density, iceDensity);
        } else {
            SLog(EError, "Unsupported snow properties calculation requested");
        }
        sigmaS = sigmaT - sigmaA;
        singleScatteringAlbedo = sigmaS / sigmaT;
    }
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
