/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/consttexture.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/core/random.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/*!
The Lambertian material represents a one-sided ideal diffuse material
with the specified amount of reflectance. Optionally, a texture map may
be applied. If no extra information is provided, the material will revert to 
the default of uniform 50% reflectance.

Seen from the back side, this material will appear completely black.

\verbatim
<bsdf type="hanrahan-krueger">
    <srgb name="reflectance" value="#a4da85"/>
</bsdf>
\endverbatim
*/
class HanrahanKrueger : public BSDF {
public:
	HanrahanKrueger(const Properties &props) 
		: BSDF(props) {

        Spectrum defaultSigmaS, defaultSigmaA;

        defaultSigmaA.fromLinearRGB(0.0014f, 0.0025f, 0.0142f);
        defaultSigmaS.fromLinearRGB(0.7f, 1.22f, 1.9f);

        if (props.hasProperty("densityMultiplier"))
            m_sizeMultiplier = props.getFloat("densityMultiplier");
        else
            m_sizeMultiplier = props.getFloat("sizeMultiplier", 1);
        m_sigmaA = props.getSpectrum("sigmaA", defaultSigmaA);
        m_sigmaS = props.getSpectrum("sigmaS", defaultSigmaS);
        m_sigmaA *= m_sizeMultiplier;
        m_sigmaS *= m_sizeMultiplier;
		/* Asymmetry parameter of the phase function */
		m_g = props.getFloat("g", 0);
        /* Relative index of refraction */
		m_eta = props.getFloat("eta", 1.32);

		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;
	}

	HanrahanKrueger(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;

        configure();
	}

	virtual ~HanrahanKrueger() {
		delete[] m_type;
	}

    void configure() {
        /* we need the scene for intersection tests */
        //m_scene = static_cast<Scene *>(getResource("scene"));
        /* Calculate extinction coefficient */
        m_sigmaT = m_sigmaA + m_sigmaS;
        /* get the longest mean free path */
        m_invSigmaTMin = 1.0f / m_sigmaT.min();
        m_invSigmaT = m_sigmaT.pow(-1.0f);
        m_negSigmaT = m_sigmaT * (-1.0f);

        /* Calculate albedo */
        m_albedo = (m_sigmaS/m_sigmaT).max();
    }

    /**
     *  Evaluate the Henyey-Greenstein phase function for two vectors with
        an asymmetry value g.  v1 and v2 should be normalized and g should 
     *  be in the range (-1, 1).  Negative values of g correspond to more
     *  back-scattering and positive values correspond to more forward scattering.
     */
    Float hgPhaseFunction(const Vector& v1, const Vector& v2, Float g) const {
	    Float costheta = dot(-v1, v2);
        Float gSq = g*g;
        Float num = 1.0 - gSq;
        Float den = std::pow(1.0 + gSq - 2.0 * g * costheta, 1.5);

    	return num / den;
    }

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		//return m_reflectance->getValue(its);
        return Spectrum(0.0f);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return Spectrum(0.0f);

        return radiance(bRec) * INV_PI;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

    /**
     * Computes the single-scattering radiance.
     */
    Spectrum radiance(const BSDFQueryRecord &bRec) const {
        /* cosines of input and output directions */
        const Float cos_wi = Frame::cosTheta(bRec.wi);
        const Float cos_wo = Frame::cosTheta(bRec.wo);
        const Float cos_wo_abs = std::abs(cos_wo);

        Float oneovereta = 1.0f / m_eta;
        Float oneoveretaSq = oneovereta * oneovereta;;

        /* Using Snell's law, calculate the squared sine of the
         * angle between the normal and the transmitted ray */
        Float sinTheta2Sqr = oneoveretaSq * Frame::sinTheta2(bRec.wi);

        if (sinTheta2Sqr > 1.0f) /* Total internal reflection! */
            return Spectrum(1.0f);

        /* Compute the cosine, but guard against numerical imprecision */
        Float cosTheta2 = std::sqrt(std::max((Float) 0.0f, 1.0f - sinTheta2Sqr));
        /* With cos(N, transmittedRay) on tap, calculating the 
         * transmission direction is straightforward. */
        Vector localTo = Vector(-oneovereta*bRec.wi.x, -oneovereta*bRec.wi.y, -cosTheta2);
        Vector to = normalize( bRec.its.toWorld( localTo ));

        /* importance sampling norminator */
        Random* random = m_random.get();
        if (random == NULL) {
            random = new Random();
            m_random.set(random);
        }
        Float sample = random->nextFloat();
        if (sample < 0.001)
            sample = 0.001;
        const Float ran = - std::log( sample );
        /* so' with max. maen frea path */
        const Float soPrimeMin =  m_invSigmaTMin * ran;

        /* Get sample point on refracted ray in world coordinates */
        const Point &xi = bRec.its.p;
        const Point3 xsamp = xi + to * soPrimeMin;

        /* Calculate siPrime and soPrime */

        /* Indireclty find intersection of light with surface xo by using the
         * triangle xi, xo, xamp with angles ai, ao, asamp. By using the
         * height/z-difference between xi and xsamp, we can calculate
         * si = h/(sin ao). si is the distance from sample point in surface
         * to light entering point. If gamma is the angle between normal an wo,
         * then ao = 90 degree - gamma. The sine of ao eqals the sine of
         * (90 - gamma), which again is (sin 90 * cos gamma - cos 90 * sin gamma)
         * This can be reduced to cos gamma.
         */
        const Float si = std::abs(xi.z - xsamp.z) / cos_wo;

        /* so' over whole spectrum */
        const Float term = 1.0f - (cos_wo_abs * cos_wo_abs);
        const Float siPrime = si * cos_wo_abs / sqrt(1.0f - oneoveretaSq * term);

        /* Calculate combined transmission coefficient */
        const Float G = std::abs(cosTheta2) / cos_wo_abs;
        const Spectrum sigmaTc = m_sigmaT + m_sigmaT * G;

        /* Calculate Fresnel trensmission T = 1- R */
        const Float Ft1 = 1.0f - fresnel(cos_wo, 1.0f, m_eta);
        const Float Ft2 = 1.0f - fresnel(cos_wi, 1.0f, m_eta);
        const Float F = Ft1 * Ft2;

        /* Query phase function */
        const Float p = hgPhaseFunction(bRec.wi, bRec.wo, m_g);

        const Spectrum siTerm = (m_negSigmaT * siPrime).exp();
        /* Actually the soTerm would be e^(-sPrime_o * sigmaT), but
         * this could be reduced to e^(ran) since sPrime_o = -ran/sigmaT. */
        const Spectrum soTerm = Spectrum( exp(ran) );

        Spectrum Lo = (m_sigmaS * F * p / sigmaTc) * siTerm * soTerm;
        return Lo;
    }

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

        return radiance(bRec);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
        pdf = Frame::cosTheta(bRec.wo) * INV_PI;

        return radiance(bRec) * INV_PI;
	}
		
	void addChild(const std::string &name, ConfigurableObject *child) {
//		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "reflectance") {
//			m_reflectance = static_cast<Texture *>(child);
//			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
//		} else {
			BSDF::addChild(name, child);
//		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
        stream->writeFloat(m_sizeMultiplier);
        m_sigmaA.serialize(stream);
        m_sigmaS.serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Hanrahan-Krueger[" << std::endl
            << "  sigmaA = " << m_sigmaA.toString() << "," << std::endl
            << "  sigmaS = " << m_sigmaS.toString() << "," << std::endl
            << "  sigmaT = " << m_sigmaT.toString() << "," << std::endl
            << "]";

		return oss.str();
	}

//	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
    Spectrum m_sigmaS;
    Spectrum m_sigmaA;
    Spectrum m_sigmaT;
    Spectrum m_invSigmaT;
    Spectrum m_negSigmaT;
    Float m_invSigmaTMin; 
    Float m_sizeMultiplier;
    Float m_g;
    Float m_eta;
    Float m_albedo;
    mutable ThreadLocal<Random> m_random;
    ref<Scene> m_scene;
};

// ================ Hardware shader implementation ================ 
//
//class HanrahanKruegerShader : public Shader {
//public:
//	HanrahanKruegerShader(Renderer *renderer, const Texture *reflectance) 
//		: Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
//		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
//	}
//
//	bool isComplete() const {
//		return m_reflectanceShader.get() != NULL;
//	}
//
//	void cleanup(Renderer *renderer) {
//		renderer->unregisterShaderForResource(m_reflectance.get());
//	}
//
//	void putDependencies(std::vector<Shader *> &deps) {
//		deps.push_back(m_reflectanceShader.get());
//	}
//
//	void generateCode(std::ostringstream &oss,
//			const std::string &evalName,
//			const std::vector<std::string> &depNames) const {
//		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
//			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
//			<< "    	return vec3(0.0);" << endl
//			<< "    return " << depNames[0] << "(uv) * 0.31831;" << endl
//			<< "}" << endl
//			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
//			<< "    return " << evalName << "(uv, wi, wo);" << endl
//			<< "}" << endl;
//	}
//
//	MTS_DECLARE_CLASS()
//private:
//	ref<const Texture> m_reflectance;
//	ref<Shader> m_reflectanceShader;
//};
//
//Shader *HanrahanKrueger::createShader(Renderer *renderer) const { 
//	return new HanrahanKruegerShader(renderer, m_reflectance.get());
//}
//
//MTS_IMPLEMENT_CLASS(HanrahanKruegerShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(HanrahanKrueger, false, BSDF)
MTS_EXPORT_PLUGIN(HanrahanKrueger, "Hanrahan-Krueger BRDF")
MTS_NAMESPACE_END
