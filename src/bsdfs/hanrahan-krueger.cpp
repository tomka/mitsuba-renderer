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
class HanrahanKruegar : public BSDF {
public:
	HanrahanKruegar(const Properties &props) 
		: BSDF(props) {

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

        configure();
	}

	HanrahanKruegar(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;

        configure();
	}

	virtual ~HanrahanKruegar() {
		delete[] m_type;
	}

    void configure() {
        /* Calculate extinction coefficient */
        m_sigmaT = m_sigmaA + m_sigmaS;
        /* Calculate albedo */
        m_albedo = (m_sigmaS/m_sigmaT).max();

        m_random = new Random();
    }

    /**
     *  Evaluate the Henyey-Greenstein phase function for two vectors with
        an asymmetry value g.  v1 and v2 should be normalized and g should 
     *  be in the range (-1, 1).  Negative values of g correspond to more
     *  back-scattering and positive values correspond to more forward scattering.
     */
    Float hgPhaseFunction(Vector v1, Vector v2, Float g) {
	    Float costheta = dot(-v1, v2);
        Float num = 1.0 - g*g;
        Float den = std::pow(1.0 + g*g - 2.0 * g * costheta, 1.5);

    	return num / den;
    }

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_reflectance->getValue(its);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
//		if (!(bRec.typeMask & m_combinedType)
//			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
//			return Spectrum(0.0f);
//
//		return m_reflectance->getValue(bRec.its) * INV_PI;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

        Vector n = bRec.its.shFrame.n;
        Float cos_wo = dot(bRec.wo, n);
        Float cos_wi = dot(bRec.wi, n);
        /* Calculate Fresnel trensmission T = 1- R */
        Float Ft1 = 1 - fresnel(cos_wo, 1.0, m_eta);
        Float Ft2 = 1 - fresnel(cos_wi, 1.0, m_eta);
        Float F = Ft1 * Ft2;
        /* Query phase function */
        Float p = hgPhaseFunction(bRec.wo, bRec.wi, m_g);
        /* Calculate combined transmission coefficient */
        G = std::abs(cos_wo) / std::abs(cos_wi);
        Spectrum sigmaTc = sigmaT + G * sigmaT;
        /* Calculate siPrime and soPrime */
        Float soPrime = std::log( m_random.nextFloat() ) / m_sigmaT;
        Float var1 = 1 - (1 / (m_eta * m_eta));
        Float ver2 = 1 - (cos_wi * cos_wi);
        Float siPrime = si * cos_wi / sqrt(var1 * var2);

        Spectrum Lo = (m_sigmaS * F * p / sigmaTc)
            * exp(-siPrime * m_sigmaT) * exp(-soPrime * m_sigmaT);

        return Lo;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
//		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
//			return Spectrum(0.0f);
//		bRec.wo = squareToHemispherePSA(bRec.sample);
//		bRec.sampledComponent = 0;
//		bRec.sampledType = EDiffuseReflection;
//		pdf = Frame::cosTheta(bRec.wo) * INV_PI;
//		return m_reflectance->getValue(bRec.its) * INV_PI;
	}
		
	void addChild(const std::string &name, ConfigurableObject *child) {
//		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "reflectance") {
//			m_reflectance = static_cast<Texture *>(child);
//			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
//		} else {
//			BSDF::addChild(name, child);
//		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
        stream->writeFloat(m_sizeMultiplier);
        m_sigmaA.serialize(stream);
        m_sigmaS.serialize(stream);
        manager->serialize(stream, m_phaseFunction.get());
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
    Float m_sizeMultiplier;
    Float m_g;
    Float m_eta;
    Float m_albedo;
    Random m_random;
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
