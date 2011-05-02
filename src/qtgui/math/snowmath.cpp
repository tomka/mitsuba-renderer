#ifndef SNOWMATH_H
#define SNOWMATH_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

Spectrum getSigmaAofIce() {
    InterpolatedSpectrum smoothSigmaA(7);
    smoothSigmaA.appendSample(400, 0.085);
    smoothSigmaA.appendSample(450, 0.043);
    smoothSigmaA.appendSample(500, 0.048);
    smoothSigmaA.appendSample(550, 0.071);
    smoothSigmaA.appendSample(600, 0.120);
    smoothSigmaA.appendSample(650, 0.276);
    smoothSigmaA.appendSample(700, 0.520);

    Spectrum sigmaA;
    sigmaA.fromSmoothSpectrum(&smoothSigmaA);

    return sigmaA;
}

inline Spectrum getDiffusionZr(Spectrum &sigmaTPrime) {
    return 1 / sigmaTPrime;
}

inline Spectrum getDiffusionZv(Spectrum &zr, Float A) {
    return zr * (1 + ((4*A) / 3));
}

inline Float getA(Float Fdr) {
    return (1 + Fdr) / (1 - Fdr);
}

inline Float getFdr(Float n) {
        if (n > 1)
            return (-1.4399 / (n*n)) + (0.7099 / n) + 0.6681 + 0.0636 * n;
        else
            return -0.4399 + (0.7099 / n) - (0.3319 / (n*n)) + (0.0636 / (n*n*n));
}

inline Spectrum getAlbedo(Spectrum &sigmaAIce, Float d) {
        Spectrum sigmaAIceSq(sigmaAIce);
        sigmaAIceSq.sqrt();
        return Spectrum(1.0f) - (5.96 * sigmaAIceSq * sqrt(d));
}

inline Spectrum getSingleScatteringAlbedo(Spectrum &sigmaA_ice, Float d) {
        return Spectrum(1.0f) - (0.84 * sigmaA_ice * d);
}

inline Float getNumberDensity(Float d, Float rho, Float rhoIce) {
        return (6.0f / (M_PI * d * d * d)) * (rho / rhoIce);
}

inline Spectrum getSigmaA(Spectrum &absIce, Float rho, Float rhoIce) {
        return 1.26f * absIce * (rho / rhoIce);
}

inline Spectrum getSigmaT(Float d, Float rho, Float rhoIce) {
        // Calculate geometrical cross-section
        Float G = M_PI * d * d * 0.25f;
        // Calculate extinction cross-section
        Float Cext = 2.0f * G;
        // Get number density
        Float N = getNumberDensity(d, rho, rhoIce);
        return Spectrum(Cext * N);
}

inline Spectrum getAsymptoticExtCoeff(Spectrum &absCoeffIce, Float d, Float rho, Float rhoIce) {
        Spectrum absCoeffIceSq(absCoeffIce);
        absCoeffIceSq.sqrt();
        return 0.845f * absCoeffIceSq * (1 / sqrt(d)) * (rho/rhoIce);
}

inline Spectrum getBarkstromExtCoeff(Spectrum &absCoeffIce, Float d, Float rho, Float rhoIce, Float v0 = 5.80) {
    Spectrum Ks = getAsymptoticExtCoeff(absCoeffIce, d, rho, rhoIce);
    return v0 * Ks;
}

inline Spectrum getBarkstromAbsCoeff(Spectrum &singleScatAlbedo, Spectrum &extCoeff) {
    return (Spectrum(1.0f) - singleScatAlbedo) * extCoeff;
}

inline Spectrum getReducedScatterCoeff(Spectrum &sigmaS, Float g) {
    return sigmaS * (1 - g);
}

MTS_NAMESPACE_END

#endif // SNOWMATH_H
