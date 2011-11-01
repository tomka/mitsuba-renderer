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

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gpuprogram.h>

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
	}

	Wiscombe(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager), one(1.0f) {
        m_g = stream->readFloat();
        m_depth = stream->readFloat();
        m_w0 = Spectrum(stream);
        m_sigmaT = Spectrum(stream);
	}

	virtual ~Wiscombe() {
	}

    void configure() {
		m_usesRayDifferentials = false;

        m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(Sampler::m_theClass, Properties("independent"))); 

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

		m_components.clear();
		m_components.push_back(EDiffuseReflection | EFrontSide);
		BSDF::configure();
    }

	Spectrum getDiffuseReflectance(const Intersection &its) const {
	    //return reflectance(mu0, muPrime) * mu0;
        return Spectrum(0.0f);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return Spectrum(0.0f);

        Float mu0 = Frame::cosTheta(bRec.wo);
        Float muPrime = Frame::cosTheta(bRec.wi);
		return reflectance(mu0, muPrime) * INV_PI;
	}

    /**
     * Calculate the reflectane for a specific value for u0.
     * With u0 = cos( incident beam ).
     */
    Spectrum reflectance(Float mu0, Float muPrime) const {
        Float b = 1.07 * mu0 - 0.84;
        Float fBar = (3 / (3 - b)) * (1 + b * (muPrime - 1));
        Spectrum R = albedo(mu0) * fBar * INV_PI;
        return R;
    }

    /**
     * Calculate the albedo for a specific value for u0.
     * With u0 = cos( incident beam ).
     */
    Spectrum albedo(Float mu0) const {
        return ( (m_wStar) / (one + m_P) ) * ( (one - m_xi * mu0 * m_bStar) / (one + m_xi * mu0) );
    }

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return 0.0f;

		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

        Float mu0 = Frame::cosTheta(bRec.wo);
        Float muPrime = std::abs(Frame::cosTheta(bRec.wi));
		return reflectance(mu0, muPrime) / mu0;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = Frame::cosTheta(bRec.wo) * INV_PI;

        Float mu0 = Frame::cosTheta(bRec.wo);
        Float muPrime = std::abs(Frame::cosTheta(bRec.wi));
		return reflectance(mu0, muPrime);
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

	Shader *createShader(Renderer *renderer) const;

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
    // a sampler
    ref<Sampler> m_sampler;
};

// ================ Hardware shader implementation ================ 

class WiscombeShader : public Shader {
public:
	WiscombeShader(Renderer *renderer,
        Spectrum wStar, Spectrum bStar, Spectrum P, Spectrum xi) 
		: Shader(renderer, EBSDFShader),
          m_wStar(wStar), m_bStar(bStar), m_P(P), m_xi(xi) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss	<< "uniform vec3 " << evalName << "_wStar;" << endl
			<< "uniform vec3 " << evalName << "_bStar;" << endl
			<< "uniform vec3 " << evalName << "_P;" << endl
			<< "uniform vec3 " << evalName << "_xi;" << endl
            << endl
            << "const float " << evalName << "_PI = 3.14159265358979323846264;" << endl
            << "const float " << evalName << "_invpi = 1.0 / " << evalName << "_PI;" << endl
            << "const vec3  " << evalName << "_one = vec3(1.0,1.0,1.0);" << endl 
            << endl
            << "vec3 " << evalName << "_albedo(float mu0) {" << endl
            << "    vec3 albedo = (" << evalName << "_wStar / (vec3(1.0) + " << evalName << "_P))" << endl
            << "       * ( (vec3(1.0) - " << evalName << "_xi * mu0 * " << evalName << "_bStar) / (vec3(1.0) + "
            << "      " << evalName << "_xi * mu0) );" << endl
            << "    return albedo;" << endl
            << "}" << endl
            << endl
		    << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
            << "    float mu0 = dot(wi, normal);" << endl
            << "    float muPrime = dot(wo, normal);" << endl
            << "    float b = 1.07 * mu0 - 0.84;" << endl
            << "    float fBar = (3 / (3 -b)) * (1 + b * (muPrime - 1));" << endl
            << "    vec3 albedo = " << evalName << "_albedo(mu0);" << endl
            << "    return albedo * fBar * " << evalName << "_invpi;" << endl
			<< "}" << endl
            << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_wStar", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_bStar", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_P", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_xi", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_wStar);
		program->setParameter(parameterIDs[1], m_bStar);
		program->setParameter(parameterIDs[2], m_P);
		program->setParameter(parameterIDs[3], m_xi);
	}
	MTS_DECLARE_CLASS()
private:
    Spectrum m_wStar, m_bStar, m_P, m_xi;
};

Shader *Wiscombe::createShader(Renderer *renderer) const { 
	return new WiscombeShader(renderer, m_wStar, m_bStar, m_P, m_xi);
}

MTS_IMPLEMENT_CLASS(WiscombeShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Wiscombe, false, BSDF)
MTS_EXPORT_PLUGIN(Wiscombe, "Wiscombe snow BRDF")
MTS_NAMESPACE_END
