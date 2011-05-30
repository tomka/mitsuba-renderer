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

#if !defined(__DIRECT_HW_H)
#define __DIRECT_HW_H

#include <mitsuba/hw/renderer.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * Manage the direct rendering of a scene by making use of
 * GPU programs. Real-Rime subsurface scattering is supported.
 */
class MTS_EXPORT_HW DirectShaderManager : public Object {
public:
	DirectShaderManager(const Scene *scene, Renderer *renderer);

	/// To be called once before use
	void init();

	/// Prepare for rendering a material with BSDF 'bsdf'
	void configure(const BSDF *bsdf, const Luminaire *luminaire,
            const Point &camPos, bool faceNormals);

	/// Draw the background if there is an environment luminaire
	void drawBackground(const Transform &clipToWorld, const Point &camPos);

	/// Release bound resources
	void unbind();

	/// Return all bound triangle meshes
	inline const std::vector<const TriMesh *> &getMeshes() const { return m_meshes; }

	/// Return the shadow cube map for debugging purposes
	inline GPUTexture *getShadowMap() { return m_shadowMap; }

	/// Should the shadow map generation be done in a single pass? (requires geometry shader support)
	inline void setSinglePass(bool singlePass) { m_singlePass = singlePass; }
	
	/// Is shadow map generation generation performed in a single pass?
	inline bool isSinglePass() const { return m_singlePass; }

	/// Set whether or not surfaces are drawn assumed to be diffuse
	inline void setDiffuseReceivers(bool diffuseReceivers) { m_diffuseReceivers = diffuseReceivers; }

	/// Return whether or not surfaces are assumed to be diffuse
	inline bool getDiffuseReceivers() const { return m_diffuseReceivers; }

	/// Set the current shadow map resolution
	inline void setShadowMapResolution(int resolution) { m_shadowMapResolution = resolution; }

	/// Return the current shadow map resolution
	inline int getShadowMapResolution() const { return m_shadowMapResolution; }

	/// Set the max. shadow map far plane distance
	inline void setMaxClipDist(Float maxClipDist) { m_maxClipDist = maxClipDist; }
	
	/// Return the max. shadow map far plane distance
	inline Float getMaxClipDist() const { return m_maxClipDist; }

	/// Set the clamping distance
	inline void setClamping(Float clamping) { m_clamping = clamping; }

	/// Return the clamping distance
	inline Float getClamping() const { return m_clamping; }

	/// Return the associated scene
	inline const Scene *getScene() const { return m_scene.get(); }

	/// To be called once before destruction
	void cleanup();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~DirectShaderManager();
private:
    struct DirectDependencyNode {
        Shader *shader;
        std::vector<DirectDependencyNode> children;
        std::vector<int> parameterIDs;

        inline DirectDependencyNode(const DirectDependencyNode &node)
            : shader(node.shader), children(node.children), parameterIDs(node.parameterIDs) {
        }

        inline DirectDependencyNode(Shader *shader = NULL) : shader(shader) {
            if (shader == NULL)
                return;
            std::vector<Shader *> deps;
            shader->putDependencies(deps);
            for (std::vector<Shader *>::iterator it = deps.begin();
                it != deps.end(); ++it)
                children.push_back(DirectDependencyNode(*it));
        }

        std::string recursiveGenerateCode(std::ostringstream &oss, int &id) const {
            std::vector<std::string> depNames;
            for (size_t i=0; i<children.size(); ++i)
                depNames.push_back(children[i].recursiveGenerateCode(oss, id));
            std::string evalName = formatString("shader_%i", id++);
            shader->generateCode(oss, evalName, depNames);
            oss << endl;
            return evalName;
        }

        void recursiveResolve(GPUProgram *program, int &id) {
            std::vector<std::string> depNames;
            for (size_t i=0; i<children.size(); ++i)
                children[i].recursiveResolve(program, id);

            std::string evalName = formatString("shader_%i", id++);
            shader->resolve(program, evalName, parameterIDs);
        }

        void recursiveBind(GPUProgram *program, const DirectDependencyNode &targetNode, int &textureUnitOffset) {
            for (size_t i=0; i<children.size(); ++i)
                children[i].recursiveBind(program, targetNode.children[i], textureUnitOffset);
            shader->bind(program, targetNode.parameterIDs, textureUnitOffset);
        }

        void recursiveUnbind() {
            shader->unbind();
            for (size_t i=0; i<children.size(); ++i)
                children[i].recursiveUnbind();
        }

        inline void toString(std::ostringstream &oss) const {
            oss << shader->getClass()->getName();
            if (children.size() > 0) {
                oss << '{';
                for (size_t i=0; i<children.size(); ++i) {
                    children[i].toString(oss);
                    if (i+1<children.size())
                        oss << ',';
                }
                oss << "}";
            }
        }
    };

    struct DirectProgramConfiguration {
        int param_lightPos, param_lightDir, param_lightColor;
        int param_lightAperture;

        inline void toString(std::ostringstream &oss) const {
            oss << "vpl=";
            //vpl.toString(oss);
            oss << ", bsdf=";
            //bsdf.toString(oss);
            //if (hasLuminaire) {
            //    oss << ", luminaire=";
            //    luminaire.toString(oss);
            //}
            //if (faceNormals)
            //    oss << ", faceNormals";
        }
    };

    struct ProgramAndConfiguration {
        GPUProgram *program;
        DirectProgramConfiguration config;

        inline ProgramAndConfiguration() : program(NULL) {
        }

        inline ProgramAndConfiguration(GPUProgram *program, 
            const DirectProgramConfiguration &config)
            : program(program), config(config) {
        }
    };

	/* General */
	ref<const Scene> m_scene;
	ref<Renderer> m_renderer;
	Float m_clamping, m_minDist;
	Float m_maxClipDist;
	bool m_initialized;
	
	/* Shadow mapping related */
	ref<GPUProgram> m_shadowProgram;
	ref<GPUProgram> m_altShadowProgram;
	int m_shadowProgramParam_cubeMapTransform[6];
	int m_shadowProgramParam_depthVec[6];
	int m_altShadowProgramParam_cubeMapTransform;
	int m_altShadowProgramParam_depthVec;
	ref<GPUTexture> m_shadowMap;
	Float m_nearClip, m_invClipRange;
	int m_shadowMapResolution;
	bool m_singlePass;
	bool m_diffuseReceivers;

	/* Rendering related */
    std::map<std::string, ProgramAndConfiguration> m_programs;
    ProgramAndConfiguration m_current;
    DirectProgramConfiguration m_targetConfig;
	ref<GPUProgram> m_backgroundProgram;
    DirectDependencyNode m_backgroundDependencies; 
	std::vector<const TriMesh *> m_meshes;

public:
    ref<GPUProgram> m_lightViewProgram;
    int param_lightPos, param_lightDir, param_lightColor,
        param_lightAperture, param_lightAlbedoTex;
    
    ref<GPUProgram> m_cameraViewProgram;
    ref<GPUProgram> m_expandSilhouetteProgram;
    int param_expandViewTex;

    ref<GPUProgram> m_renderSplatsProgram;
    int attrib_renderSplatsBillboardOffset, param_renderSplatsViewSurfacePos,
        param_renderSplatsTranslucencyTex, param_renderSplatsBillboardRadius;

    ref<GPUProgram> m_finalContributionProgram;
    int param_finalContribSubSurf, param_finalContribAlbedoTex,
        param_finalContribSampleScale, param_finalContribLightAperture,
        param_finalContribLightSpecColor, param_finalContribLightDir,
        param_finalContribLightPos;
};



MTS_NAMESPACE_END

#endif /* __DIRECT_HW_H */
