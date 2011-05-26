#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/direct.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

DirectShaderManager::DirectShaderManager(const Scene *scene, Renderer *renderer)
     : m_scene(scene), m_renderer(renderer), m_clamping(0.1f),
       m_maxClipDist(std::numeric_limits<Float>::infinity()), m_initialized(false), 
       m_shadowMapResolution(512), m_singlePass(false), 
       m_diffuseReceivers(false) {
}

DirectShaderManager::~DirectShaderManager() {
    if (m_initialized)
        cleanup();
}


void DirectShaderManager::init() {
    if (m_renderer->getCapabilities()->isSupported(RendererCapabilities::EGeometryShaders)) {
        m_shadowProgram = m_renderer->createGPUProgram("Shadow Program");

        /* In the shadow vertex program we want to stay in object space */
        m_shadowProgram->setSource(GPUProgram::EVertexProgram,
            "void main() {\n"
            "   gl_Position = gl_Vertex;\n"
            "}"
        );

        m_shadowProgram->setSource(GPUProgram::EGeometryProgram,
            "#version 120\n"
            "#extension GL_EXT_geometry_shader4 : enable\n"
            "\n"
            "uniform mat4 cubeMapTransform[6];\n"
            "uniform vec4 depthVec[6];\n"
            "varying float depth;\n"
            "\n"
            "void main() {\n"
            "   depth = 0;\n" // avoid a (incorrect?) warning
            "   for (int side = 0; side < 6; side++) {\n"
            "       gl_Layer = side;\n"
            "       for (int i = 0; i < gl_VerticesIn; i++) {\n"
            "           gl_Position = cubeMapTransform[side] * gl_PositionIn[i];\n"
            "           depth = dot(depthVec[side], gl_PositionIn[i]);\n"
            "           EmitVertex();\n"
            "       }\n"
            "       EndPrimitive();\n"
            "   }\n"
            "}\n"
        );

        m_shadowProgram->setSource(GPUProgram::EFragmentProgram,
            "#version 120\n"
            "varying float depth;\n"
            "void main() {\n"
            "   float dx = dFdx(depth), dy = dFdy(depth);"
            "   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
            "}\n"
        );

        /* Six output triangles per input triangle */
        m_shadowProgram->setMaxVertices(18);
        m_shadowProgram->init(); 

        for (int i=0; i<6; ++i) {
            m_shadowProgramParam_cubeMapTransform[i] =
                m_shadowProgram->getParameterID(formatString("cubeMapTransform[%i]", i));
            m_shadowProgramParam_depthVec[i] =
                m_shadowProgram->getParameterID(formatString("depthVec[%i]", i));
        }
    }

    m_altShadowProgram = m_renderer->createGPUProgram("Alternative Shadow Program");
    m_altShadowProgram->setSource(GPUProgram::EVertexProgram,
        "uniform mat4 cubeMapTransform;\n"
        "uniform vec4 depthVec;\n"
        "varying float depth;\n"
        "void main() {\n"
        "   gl_Position = cubeMapTransform * gl_Vertex;\n"
        "   depth = dot(depthVec, gl_Vertex);\n"
        "}\n"
    );

    m_altShadowProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "varying float depth;\n"
        "void main() {\n"
        "   float dx = dFdx(depth), dy = dFdy(depth);"
        "   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
        "}\n"
    );

    m_altShadowProgram->init();

    m_altShadowProgramParam_cubeMapTransform =
        m_altShadowProgram->getParameterID("cubeMapTransform");
    m_altShadowProgramParam_depthVec =
        m_altShadowProgram->getParameterID("depthVec");

    const std::vector<Shape *> shapes = m_scene->getShapes();
    const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();

    /* Try to get shader implementations for the shapes BSDF to
     * register them with the renderer.  */
    for (size_t i=0; i<shapes.size(); ++i) {
        ref<TriMesh> triMesh = shapes[i]->createTriMesh();
        if (!triMesh)
            continue;
        m_renderer->registerGeometry(triMesh);
        Shader *shader = m_renderer->registerShaderForResource(triMesh->getBSDF());
        if (shader != NULL && !shader->isComplete())
            m_renderer->unregisterShaderForResource(triMesh->getBSDF());
        triMesh->incRef();
        m_meshes.push_back(triMesh);
    }

    /* Register shader implementations for luminaires with renderer. */
    for (size_t i=0; i<luminaires.size(); ++i)
        m_renderer->registerShaderForResource(luminaires[i]);

    if (m_scene->hasBackgroundLuminaire() &&
        m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire()) != NULL) {
        Shader *shader = m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire());
        m_backgroundDependencies = DirectDependencyNode(shader);
        int id = 0;
        std::ostringstream oss;
        std::string evalName = m_backgroundDependencies.recursiveGenerateCode(oss, id);

        m_backgroundProgram = m_renderer->createGPUProgram("Background program");
        m_backgroundProgram->setSource(GPUProgram::EVertexProgram,
            "uniform mat4 clipToWorld;\n"
            "varying vec3 d;\n"
            "void main() {\n"
            "   gl_Position = ftransform();\n"
            "   vec4 tmp = clipToWorld * (gl_ModelViewProjectionMatrix * gl_Vertex);\n"
            "   d = tmp.xyz/tmp.w;"
            "}\n"
        );

        oss << "varying vec3 d;" << endl
            << "uniform vec3 camPos;" << endl
            << "void main() {" << endl
            << "  gl_FragColor.rgb = " << evalName << "_background(normalize(d - camPos));" << endl
            << "  gl_FragColor.a = 1.0;" << endl
            << "}" << endl;

        m_backgroundProgram->setSource(GPUProgram::EFragmentProgram, oss.str());
        m_backgroundProgram->init();

        id = 0;
        m_backgroundDependencies.recursiveResolve(m_backgroundProgram, id);
    }

    m_initialized = true;
}

void DirectShaderManager::drawBackground(const Transform &clipToWorld, const Point &camPos) {
    if (m_backgroundProgram == NULL)
        return;
    int textureUnitOffset = 0;  
    m_backgroundProgram->bind();
    m_backgroundDependencies.recursiveBind(m_backgroundProgram, 
        m_backgroundDependencies, textureUnitOffset);
    m_backgroundProgram->setParameter("clipToWorld", clipToWorld, false);
    m_backgroundProgram->setParameter("camPos", camPos, false);
    m_renderer->blitQuad(false);
    m_backgroundProgram->unbind();
    m_backgroundDependencies.recursiveUnbind();
}

void DirectShaderManager::unbind() {
    /*
    if (m_current.program && m_current.program->isBound()) {
        m_targetConfig.unbind();
        m_current.program->unbind();
        m_shadowMap->unbind();
    }
    */
}

void DirectShaderManager::cleanup() {
    if (m_shadowMap)
        m_shadowMap->cleanup();

    if (m_backgroundProgram) {
        m_backgroundProgram->cleanup();
        m_backgroundProgram = NULL;
    }

    const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();

    for (size_t i=0; i<m_meshes.size(); ++i) {
        m_renderer->unregisterGeometry(m_meshes[i]);
        m_renderer->unregisterShaderForResource(m_meshes[i]->getBSDF());
        m_meshes[i]->decRef();
    }
    m_meshes.clear();
    for (size_t i=0; i<luminaires.size(); ++i)
        m_renderer->unregisterShaderForResource(luminaires[i]);
    m_initialized = false;
}

MTS_IMPLEMENT_CLASS(DirectShaderManager, false, Object)
MTS_NAMESPACE_END

