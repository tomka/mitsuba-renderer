#ifndef SNOWMATH_H
#define SNOWMATH_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 *  Create a spectrum of the absorbtion coefficient of ice, sigmaA.
 *  Lambda is given in nm and sigmaA in 1/m.
 */
Spectrum getSigmaAofIce();

/**
 * Calculate the distance of a virtual point light of the dipole
 * approximation below the surface. The reduced extinction
 * coefficient sigmaTPrime is needed for this.
 */
inline Spectrum getDiffusionZr(const Spectrum &sigmaTPrime);

/**
 * Calculate the distance of a virtual point light of the dipole
 * approximation above the surface. The position of the counter
 * part below the surface zr is needed, as well as the change of
 * scalar irradiance (due to internal reflection at the surface) A.
 */
inline Spectrum getDiffusionZv(const Spectrum &zr, Float A);

/**      
 * Calculate the change in scalar irradiance (fluence) due to
 * internal reflection at the surface. The diffuse Fresnel term
 * Fdr is needed for this.
 */
inline Float getA(Float Fdr);

/**
 * Calculate the diffuse Fresnel term, based on the relative
 * index of refraction n. This is only needed if the indices of
 * refraction mismatch. If so, n is the ratio n1/n2.
 */
inline Float getFdr(Float n);

/**
 * Calculates the albodo of snow, based on the absorption coefficient
 * of ice and tha grain diameter d. Make sure both are basde on the same unit.
 */
Spectrum getAlbedo(const Spectrum &sigmaAIce, Float d);

/**
 * Calculates the single scattering albodo of snow, based on the absorption
 * coefficient of ice and tha grain radis d. Make sure both are basde on the
 * same unit. From: Bohren1983
 */
Spectrum getSingleScatteringAlbedo(const Spectrum &sigmaA_ice, Float d);

/**
 * Get number density of snow.
 */
Float getNumberDensity(Float d, Float rho, Float rhoIce);

/**
 * Calculates the absorption coefficient of snow. Make sure that
 * the absorption coefficient of ice 'absIce' and the density 'rho'
 * have the same unit base as well as the density of ice 'rhoIce'
 */
Spectrum getSigmaA(const Spectrum &absIce, Float rho, Float rhoIce);

/**
 * Calculates the extinction coefficient of a large particle. This is
 * basically  2 * pi * d^2 / 4 where d is the diameter.
 */
Spectrum getLargeParticleExtCoeff(Float d, Float rho, Float rhoIce);

/**
 * Calculates the asymtopic extinction coefficient of snow. Make
 * sure that the absorption coefficient of ice 'absCoeffIce', the
 * density 'rho' and the grain size 'd' have the same unit base
 * as well as the density of ice 'rhoIce'.
 */
Spectrum getAsymptoticExtCoeff(const Spectrum &absCoeffIce, Float d, Float rho, Float rhoIce);

/**
 * Calculates the real extinction coefficient as stated in Bakstrom1972.
 * It takes the same parameters as the asymptotic ext. coeff, but
 * additionally takes an eigenvalue v0. This is related to single
 * scattering albedo. Examples (w - v0): 0.99 - 5.80, 0.95 - 2.63,
 * 0.90 - 1.90, 0.80 - 1.41
 */
Spectrum getBarkstromExtCoeff(const Spectrum &absCoeffIce, Float d, Float rho, Float rhoIce, const Spectrum &v0);

Spectrum getBarkstromAbsCoeff(const Spectrum &singleScatAlbedo, Spectrum &extCoeff);

/**
 * Get the eigenvector v0 of the radiative transfer equation.
 */
Spectrum geRteEigenvector(const Spectrum &singleScatAlbedo, Float g);

Spectrum getSnowPackExtCoeff(Float rho, Float c1, Float c2);

/**
 * Calculate the reduced scattering coefficient
 */
Spectrum getReducedScatterCoeff(const Spectrum &sigmaS, Float g);

MTS_NAMESPACE_END

#endif // SNOWMATH_H
