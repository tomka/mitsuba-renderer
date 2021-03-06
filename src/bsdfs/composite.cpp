/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Composite material, represents a linear combination of 
 * one or more BRDFs.
 */
class Composite : public BSDF {
public:
	Composite(const Properties &props) 
		: BSDF(props), m_bsdfWeights(NULL), m_componentIndex(NULL), m_bsdfOffset(NULL) {
		/* Parse the weight parameter */
		std::vector<std::string> weights = 
			tokenize(props.getString("weights", ""), " ,;");
		m_bsdfCount = weights.size();
		m_bsdfWeights = new Float[m_bsdfCount];

		Float totalWeight = 0;
		char *end_ptr = NULL;
		for (size_t i=0; i<weights.size(); ++i) {
			Float weight = (Float) strtod(weights[i].c_str(), &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse the BRDF weights!");
			if (weight < 0)
				SLog(EError, "Invalid BRDF weight!");
			m_bsdfWeights[i] = weight;
			totalWeight += weight;
		}

		if (totalWeight > 1)
			Log(EWarn, "Energy conservation is violated!");
	}

	Composite(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager), m_bsdfWeights(NULL), m_componentIndex(NULL), m_bsdfOffset(NULL) {
		m_bsdfCount = stream->readSize();
		m_bsdfWeights = new Float[m_bsdfCount];
		for (size_t i=0; i<m_bsdfCount; ++i) {
			m_bsdfWeights[i] = stream->readFloat();
			BSDF *bsdf = static_cast<BSDF *>(manager->getInstance(stream));
			bsdf->incRef();
			m_bsdfs.push_back(bsdf);
		}
		configure();
	}

	virtual ~Composite() {
		for (size_t i=0; i<m_bsdfs.size(); ++i)
			m_bsdfs[i]->decRef();
		if (m_type)
			delete[] m_type;
		if (m_bsdfWeights) 
			delete[] m_bsdfWeights;
		if (m_componentIndex)
			delete[] m_componentIndex;
		if (m_bsdfOffset)
			delete[] m_bsdfOffset;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		
		stream->writeSize(m_bsdfCount);
		for (size_t i=0; i<m_bsdfCount; ++i) {
			stream->writeFloat(m_bsdfWeights[i]);
			manager->serialize(stream, m_bsdfs[i]);
		}
	}

	void configure() {
		BSDF::configure();

		m_combinedType = 0;
		m_usesRayDifferentials = false;
		m_componentCount = 0;

		if (m_bsdfs.size() != m_bsdfCount)
			Log(EError, "BSDF count mismatch: " SIZE_T_FMT " bsdfs, but specified " SIZE_T_FMT " weights",
				m_bsdfs.size(), m_bsdfCount);

		for (size_t i=0; i<m_bsdfCount; ++i)
			m_componentCount += m_bsdfs[i]->getComponentCount();

		m_pdf = DiscretePDF(m_bsdfs.size());

		if (m_type)
			delete[] m_type;
		if (m_componentIndex)
			delete[] m_componentIndex;
		if (m_bsdfOffset)
			delete[] m_bsdfOffset;

		m_type = new unsigned int[m_componentCount];
		m_componentIndex = new std::pair<int, int>[m_componentCount];
		m_bsdfOffset = new int[m_bsdfCount];
		int ctr = 0, offset = 0;

		for (size_t i=0; i<m_bsdfCount; ++i) {
			const BSDF *bsdf = m_bsdfs[i];
			m_bsdfOffset[i] = offset;

			for (int j=0; j<bsdf->getComponentCount(); ++j) {
				int componentType = bsdf->getType(j);
				m_type[offset+j] = componentType;
				m_componentIndex[ctr++] = std::make_pair((int) i, j);
			}

			m_combinedType |= bsdf->getType();
			offset += bsdf->getComponentCount();
			m_usesRayDifferentials |= bsdf->usesRayDifferentials();
			m_pdf[i] = m_bsdfWeights[i];
		}
		m_pdf.build();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum result(0.0f);
		for (size_t i=0; i<m_bsdfCount; ++i)
			result+= m_bsdfs[i]->getDiffuseReflectance(its) * m_bsdfWeights[i];
		return result;
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (bRec.component == -1) {
			for (size_t i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->f(bRec) * m_bsdfWeights[i];
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			BSDFQueryRecord bRec2(bRec);
			bRec2.component = m_componentIndex[bRec.component].second;
			return m_bsdfs[idx]->f(bRec2) * m_bsdfWeights[idx];
		}

		return result;
	}

	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (bRec.component == -1) {
			for (size_t i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->fDelta(bRec) * m_bsdfWeights[i];
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			BSDFQueryRecord bRec2(bRec);
			bRec2.component = m_componentIndex[bRec.component].second;
			return m_bsdfs[idx]->fDelta(bRec2) * m_bsdfWeights[idx];
		}

		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		Float result = 0.0f;

		if (bRec.component == -1) {
			for (size_t i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->pdf(bRec) * m_pdf[i];
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			BSDFQueryRecord bRec2(bRec);
			bRec2.component = m_componentIndex[bRec.component].second;
			return m_bsdfs[idx]->pdf(bRec2);
		}

		return result;
	}
	
	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		Float result = 0.0f;

		if (bRec.component == -1) {
			for (size_t i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->pdfDelta(bRec) * m_pdf[i];
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			BSDFQueryRecord bRec2(bRec);
			bRec2.component = m_componentIndex[bRec.component].second;
			return m_bsdfs[idx]->pdfDelta(bRec2);
		}

		return result;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		Point2 sample(_sample);
		if (bRec.component == -1) {
			int entry = m_pdf.sampleReuse(sample.x);
			Spectrum result = m_bsdfs[entry]->sample(bRec, sample);

			if (result.isZero()) // sampling failed
				return result;

			bRec.sampledComponent += m_bsdfOffset[entry];

			if (bRec.sampledType & BSDF::EDelta) {
				pdf = pdfDelta(bRec);
				return fDelta(bRec);
			} else {
				pdf = Composite::pdf(bRec);
				return f(bRec);
			}
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			int tempComponent = bRec.component;
			bRec.component = m_componentIndex[bRec.component].second;
			Spectrum result = m_bsdfs[idx]->sample(bRec, pdf, sample) * m_bsdfWeights[idx];
			bRec.component = bRec.sampledComponent = tempComponent;
			return result;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);
		if (bRec.component == -1) {
			int entry = m_pdf.sampleReuse(sample.x);
			Spectrum result = m_bsdfs[entry]->sample(bRec, sample);
			
			if (result.isZero()) // sampling failed
				return result;

			bRec.sampledComponent += m_bsdfOffset[entry];

			if (bRec.sampledType & BSDF::EDelta)
				return fDelta(bRec)/pdfDelta(bRec);
			else
				return f(bRec)/pdf(bRec);
		} else {
			/* Pick out an individual component */
			int idx = m_componentIndex[bRec.component].first;
			int tempComponent = bRec.component;
			bRec.component = m_componentIndex[bRec.component].second;
			Spectrum result = m_bsdfs[idx]->sample(bRec, sample) * m_bsdfWeights[idx];
			bRec.component = bRec.sampledComponent = tempComponent;
			return result;
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			BSDF *bsdf = static_cast<BSDF *>(child);
			m_bsdfs.push_back(bsdf);
			bsdf->incRef();
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Composite[" << endl
			<< "  weights = {";
		for (size_t i=0; i<m_bsdfCount; ++i) {
			oss << " " << m_bsdfWeights[i];
			if (i + 1 < m_bsdfCount)
				oss << ",";
		}
		oss << " }," << endl
			<< "  bsdfs = {" << endl;
		for (size_t i=0; i<m_bsdfs.size(); ++i) 
			oss << "    " << indent(m_bsdfs[i]->toString(), 2) << "," << endl;
		oss << "  }" << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	size_t m_bsdfCount;
	Float *m_bsdfWeights;
	std::pair<int, int> *m_componentIndex;
	int *m_bsdfOffset;
	std::vector<BSDF *> m_bsdfs;
	DiscretePDF m_pdf;
};

// ================ Hardware shader implementation ================ 

class CompositeShader : public Shader {
public:
	CompositeShader(Renderer *renderer, const std::vector<BSDF *> &bsdfs, Float *weights) 
		: Shader(renderer, EBSDFShader), m_bsdfs(bsdfs), m_weights(weights), m_complete(false) {
		m_bsdfShader.resize(bsdfs.size());
		for (size_t i=0; i<bsdfs.size(); ++i) {
			ref<Shader> shader = renderer->registerShaderForResource(bsdfs[i]);
			if (shader) {
				shader->incRef();
				/* At least one shader has a hardware implementation */
				m_complete = true;
			}
			m_bsdfShader[i] = shader;
		}
	}

	bool isComplete() const {
		return m_complete;
	}

	void cleanup(Renderer *renderer) {
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			renderer->unregisterShaderForResource(m_bsdfs[i]);
			if (m_bsdfShader[i])
				m_bsdfShader[i]->decRef();
		}
		m_bsdfShader.clear();
	}

	void putDependencies(std::vector<Shader *> &deps) {
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (m_bsdfShader[i])
				deps.push_back(m_bsdfShader[i]);
		}
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		Assert(m_complete);
		int ctr = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (!m_bsdfShader[i])
				continue;
			oss << "uniform float " << evalName << "_weight_" << ctr++ << ";" << endl;
		}
		oss << endl;
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return ";
		ctr = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (!m_bsdfShader[i])
				continue;
			oss << endl << "      ";
			if (ctr != 0)
				oss << "+ ";
			else
				oss << "  ";
			oss << depNames[ctr] << "(uv, wi, wo) * "
							  << evalName << "_weight_" << ctr;
			ctr++;
		}
		oss << ";" << endl << "}" << endl << endl;
		oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return ";
		ctr = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (!m_bsdfShader[i])
				continue;
			oss << endl << "      ";
			if (ctr != 0)
				oss << "+ ";
			else
				oss << "  ";
			oss << depNames[ctr] << "_diffuse(uv, wi, wo) * "
							  << evalName << "_weight_" << ctr;
			ctr++;
		}
		oss << ";" << endl << "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		int ctr = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (!m_bsdfShader[i])
				continue;
			parameterIDs.push_back(
				program->getParameterID(formatString("%s_weight_%i", evalName.c_str(), ctr++)));
		}
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		int ctr = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (!m_bsdfShader[i])
				continue;
			program->setParameter(parameterIDs[ctr++], m_weights[i]);
		}
	}


	MTS_DECLARE_CLASS()
private:
	std::vector<Shader *> m_bsdfShader;
	const std::vector<BSDF *> &m_bsdfs;
	const Float *m_weights;
	bool m_complete;
};

Shader *Composite::createShader(Renderer *renderer) const { 
	return new CompositeShader(renderer, m_bsdfs, m_bsdfWeights);
}

MTS_IMPLEMENT_CLASS(CompositeShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Composite, false, BSDF)
MTS_EXPORT_PLUGIN(Composite, "Composite BRDF")
MTS_NAMESPACE_END
