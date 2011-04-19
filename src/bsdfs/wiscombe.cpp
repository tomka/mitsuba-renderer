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

MTS_NAMESPACE_BEGIN

/*!
A BRDF model which computes the 'directional hemispherical'
reflectance of snow. It has no angular distribution, but only the
elevation angle of the sun/light source wrt. to the normal.

Seen from the back side, this material will appear completely black.

\verbatim
<bsdf type="dozier">
    <float name="g" value="0.874" />
    <spectrum name="singleScatteringAlbedo" value="" />
</bsdf>
\endverbatim
*/
class Wiscombe : public BSDF {
public:
	Wiscombe(const Properties &props) 
		: BSDF(props), one(1.0f) {
        /* use default values for a 0.3mm, values are meter based */
        Spectrum defaultSigmaT;
        defaultSigmaT.fromLinearRGB(16.4967, 6.0957, 4.6547);

        m_g = props.getFloat("g", 0.874);
        m_depth = props.getFloat("depth", 1.0);
        m_w0 = props.getSpectrum("singleScatteringAlbodo", Spectrum(0.99));
        m_sigmaT = props.getSpectrum("sigmaT", defaultSigmaT);
        configure();
	}

	Wiscombe(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager), one(1.0f) {
        m_g = stream->readFloat();
        m_depth = stream->readFloat();
        m_w0 = Spectrum(stream);
        m_sigmaT = Spectrum(stream);
        configure();
	}

	virtual ~Wiscombe() {
		delete[] m_type;
	}

    void configure() {
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;
		m_usesRayDifferentials = false;

        // spectral optial depth
        m_tau = m_sigmaT * m_depth;

        // delta-Eddington transformation of actual layer
        Float gSq = m_g * m_g;
        m_wStar = ((1 - gSq) * m_w0) / (one - gSq * m_w0);
        m_gStar = m_g / (1 + m_g);
        
        Spectrum bTmp = (one - m_wStar * m_gStar);
        m_bStar = m_gStar / bTmp;
        m_tauStar = (one - m_w0 * gSq) * m_tau;
        
        m_gamma1 = (Spectrum(7.0f) - m_wStar * (4 + 3 * m_gStar)) * 0.25;
        m_gamma2 = (one - m_wStar * (4 - 3 * m_gStar)) * (-0.25);


        Spectrum xiSq = (one - m_wStar * m_gStar) * (one - m_wStar) * 3;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
            m_xi[i] = sqrt(xiSq[i]);
        }
        
        m_P = (2 * m_xi) / ( (one - m_wStar * m_gStar) * 3);
    }

	Spectrum getDiffuseReflectance(const Intersection &its) const {
        Float u0 = dot(its.wi, its.shFrame.n);
        return reflectance(u0);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return Spectrum(0.0f);

        Float u0 = dot(bRec.wo, bRec.its.shFrame.n);
        return reflectance(u0);
	}

    /**
     * Calculate the reflectane for a specific value for u0.
     * With u0 = cos( incident beam ).
     */
    Spectrum reflectance(Float u0) const {
        return ( (m_wStar) / (one + m_P) ) * ( (one - m_xi * u0 * m_bStar) / (one + m_xi * u0) );
    }

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

        Float u0 = dot(bRec.wo, bRec.its.shFrame.n);
		return reflectance(u0);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = Frame::cosTheta(bRec.wo) * INV_PI;

        Float u0 = dot(bRec.wo, bRec.its.shFrame.n);
		return reflectance(u0);
	}
		
	void addChild(const std::string &name, ConfigurableObject *child) {
	    BSDF::addChild(name, child);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

        stream->writeFloat(m_g);
        stream->writeFloat(m_depth);
        m_w0.serialize(stream);
        m_sigmaT.serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Wiscombe[" << std::endl
            << "\tdepth=" << m_depth << std::endl
            << "\tg=" << m_g << std::endl
            << "\tsigmaT=" << m_sigmaT.toString() << std::endl
            << "\tg*=" << m_gStar << std::endl
            << "\tw0=" << m_w0.toString() << std::endl
            << "\tw*=" << m_wStar.toString() << std::endl
            << "]";
		return oss.str();
	}

//	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
    // asymetrie factor
    Float m_g;
    // the depth of the layer
    Float m_depth;
    // extinction coefficient
    Spectrum m_sigmaT;
    // single-scattering albedo
    Spectrum m_w0;
    // the spectral optical depth
    Spectrum m_tau;

    // delta-Eddington transformed layer parameters
    Spectrum m_wStar;
    Float m_gStar;
    Spectrum m_bStar;
    Spectrum m_tauStar;
    // first two two-stream parameters for approx. scattering function
    Spectrum m_gamma1, m_gamma2;
    
    Spectrum m_xi, m_P;

    // needed now and then
    const Spectrum one;
};

// ================ Hardware shader implementation ================ 

//class LambertianShader : public Shader {
//public:
//	LambertianShader(Renderer *renderer, const Texture *reflectance) 
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
//Shader *Lambertian::createShader(Renderer *renderer) const { 
//	return new LambertianShader(renderer, m_reflectance.get());
//}

MTS_IMPLEMENT_CLASS_S(Wiscombe, false, BSDF)
MTS_EXPORT_PLUGIN(Wiscombe, "Wiscombe snow BRDF")
MTS_NAMESPACE_END
