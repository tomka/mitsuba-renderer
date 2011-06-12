#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/direct.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/properties.h>

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

    /* light view program - with albedo map */
    m_lightViewProgram = m_renderer->createGPUProgram("SplatSSS LightView Program");
    m_lightViewProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "uniform vec3 lightPos;\n"
        "uniform vec3 lightColor;\n"
        "uniform vec3 lightDir;\n"
        "uniform float lightAperture;\n"
        "uniform sampler2D albedoTex;\n"
        "\n"
        "varying vec3 surfPos;\n"
        "varying vec3 normal;\n"
        "varying vec3 surfToLight;\n"
        "\n"
        "void main() {\n"
           // Vertex position in camera space
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n"
           // Vertex position in clip space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n"
        "  surfPos = (gl_TextureMatrix[0] * gl_Vertex).xyz;\n"
        "  surfToLight = lightPos-surfPos.xyz;\n"
        "\n"
        "  gl_TexCoord[0] = gl_MultiTexCoord0;\n" //out
        "  normal = gl_Normal;\n"
        "  gl_Position = vMVP;\n"
        "}\n"
    );

    m_lightViewProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "uniform vec3 lightPos;\n"
        "uniform vec3 lightColor;\n"
        "uniform vec3 lightDir;\n"
        "uniform float lightAperture;\n"
        "uniform sampler2D albedoTex;\n"
        "\n"
        "varying vec3 surfPos;\n"
        "varying vec3 normal;\n"
        "varying vec3 surfToLight;\n"
        "\n"
        "void main() {\n"
        "  vec3 surfToLightNorm = normalize(surfToLight);\n"
        "  float lightIntensity = dot(normalize(normal), surfToLightNorm);\n"
           //sqrt because light is musplipied two time by albedo
        "  vec3 surfAlbedo  = sqrt(texture2D( albedoTex, gl_TexCoord[0].st ).rgb);\n"
           //spot extinction
        "  float spot = clamp( (dot(normalize(lightDir),-surfToLightNorm)-lightAperture)/(1.0-lightAperture) ,0.0,1.0);\n"
        "  gl_FragData[0] = vec4(surfPos, 0.0);\n" //splat origin
        "  gl_FragData[1] = vec4(lightColor * spot * lightIntensity * surfAlbedo, 0.0);\n"  //splat center color
        "}\n"
    );

    // upload the program
    m_lightViewProgram->init();
    // configure parameters
    param_lightPos = m_lightViewProgram->getParameterID("lightPos", false);
    param_lightColor = m_lightViewProgram->getParameterID("lightColor", false);
    param_lightDir = m_lightViewProgram->getParameterID("lightDir", false);
    param_lightAperture = m_lightViewProgram->getParameterID("lightAperture", false);
    param_lightAlbedoTex = m_lightViewProgram->getParameterID("albedoTex", false);

    /* light view program - wiscombe warren albedo */
    PluginManager *pluginManager = PluginManager::getInstance();
    Properties props("wiscombe");
    ref<BSDF> wiscombeBsdf = static_cast<BSDF *> (pluginManager->createObject(
                    BSDF::m_theClass, props));
    wiscombeBsdf->configure();
    m_wiscombeWarrenShader = wiscombeBsdf->createShader(m_renderer);
    Assert(m_wiscombeWarrenShader != NULL);

    m_lightViewWWProgram = m_renderer->createGPUProgram("SplatSSS LightView WW Program");
    m_lightViewWWProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "uniform vec3 camPos;\n"
        "uniform vec3 lightPos;\n"
        "uniform vec3 lightColor;\n"
        "uniform vec3 lightDir;\n"
        "uniform float lightAperture;\n"
        "\n"
        "varying vec3 surfPos;\n"
        "varying vec3 normal;\n"
        "varying vec3 surfToLight;\n"
        "varying vec3 surfToCam;\n"
        "\n"
        "void main() {\n"
           // Vertex position in camera space
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n"
           // Vertex position in clip space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n"
        "  surfPos = (gl_TextureMatrix[0] * gl_Vertex).xyz;\n"
        "  surfToLight = lightPos-surfPos.xyz;\n"
        "  surfToCam = camPos-surfPos.xyz;\n"
        "\n"
        "  gl_TexCoord[0] = gl_MultiTexCoord0;\n" //out
        "  normal = gl_Normal;\n"
        "  gl_Position = vMVP;\n"
        "}\n"
    );

    std::ostringstream oss;
    oss << "#version 120\n"
        << "uniform vec3 lightPos;\n"
        << "uniform vec3 lightColor;\n"
        << "uniform vec3 lightDir;\n"
        << "uniform float lightAperture;\n"
        << "\n"
        << "varying vec3 surfPos;\n"
        << "varying vec3 normal;\n"
        << "varying vec3 surfToLight;\n"
        << "varying vec3 surfToCam;\n"
        << "\n";
    // add wiscombe warren shader code;
    m_wiscombeWarrenShader->generateCode(oss, "wiscombe", std::vector<std::string>());

    oss << "\n"
        << "void main() {\n"
        << "  vec3 surfToLightNorm = normalize(surfToLight);\n"
        << "  vec3 surfToCamNorm = normalize(surfToCam);\n"
        << "  float lightIntensity = dot(normalize(normal), surfToLightNorm);\n"
           //sqrt because light is musplipied two time by albedo
        << "  vec3 surfAlbedo  = wiscombe_albedo(lightIntensity);\n"
           //spot extinction
        << "  float spot = clamp( (dot(normalize(lightDir),-surfToLightNorm)-lightAperture)/(1.0-lightAperture) ,0.0,1.0);\n"
        << "  gl_FragData[0] = vec4(surfPos, 0.0);\n" //splat origin
        << "  gl_FragData[1] = vec4(lightColor * spot * lightIntensity * surfAlbedo, 0.0);\n"  //splat center color
        << "}\n";

    m_lightViewWWProgram->setSource(GPUProgram::EFragmentProgram, oss.str());

    // upload the program
    m_lightViewWWProgram->init();
    // configure parameters
    param_lightWWCamPos = m_lightViewWWProgram->getParameterID("camPos", false);
    param_lightWWPos = m_lightViewWWProgram->getParameterID("lightPos", false);
    param_lightWWColor = m_lightViewWWProgram->getParameterID("lightColor", false);
    param_lightWWDir = m_lightViewWWProgram->getParameterID("lightDir", false);
    param_lightWWAperture = m_lightViewWWProgram->getParameterID("lightAperture", false);
    m_wiscombeWarrenShader->resolve(m_lightViewWWProgram, "wiscombe", m_wiscombeWarrenParams);

    /* camere view program */
    m_cameraViewProgram = m_renderer->createGPUProgram("SplatSSS CameraView Program");
    m_cameraViewProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "varying vec3 surfPos;\n"
        "\n"
        "void main() {\n"
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n" //vertex position in camera space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n" //vertex position in clip space
        "  surfPos = (gl_TextureMatrix[0]*gl_Vertex).xyz;\n"  //vertex * modelMatrix
        " gl_Position = vMVP;\n" // out
        "}\n"
    );

    m_cameraViewProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "varying vec3 surfPos;\n"
        "\n"
        "void main() {\n"
        "  gl_FragColor = vec4(surfPos, 1.0);\n"  //1.0 to enable silhouette expand (to avoid black silhouette artifact)
        "}\n"
    );

    // upload the program
    m_cameraViewProgram->init();


    /* silhouette expansion program */
     m_expandSilhouetteProgram = m_renderer->createGPUProgram("SplatSSS Silhouette Expansion Program");
     m_expandSilhouetteProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "//uniform sampler2D viewSubScatTex;\n"
        "\n"
        "varying vec2 offset;\n"
        "\n"
        "void main() {\n"
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n" //vertex position in camera space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n" //vertex position in clip space
        "  offset = gl_MultiTexCoord0.zw;\n"
        "\n"
        "  gl_TexCoord[0] = gl_MultiTexCoord0;\n" // out, contains tex coord in xy and tap offset in zw
        "  gl_Position = vMVP;\n"
        "}\n"
    );

    m_expandSilhouetteProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "uniform sampler2D viewSubScatTex;\n"
        "\n"
        "varying vec2 offset;\n"
        "\n"
        "void main() {\n"
        "  vec4 center = texture2D( viewSubScatTex, gl_TexCoord[0].xy );\n"
        "  vec4 tap1   = texture2D( viewSubScatTex, gl_TexCoord[0].xy + offset );\n"
        "  vec4 tap2   = texture2D( viewSubScatTex, gl_TexCoord[0].xy - offset );\n"
        "\n"
        "  vec3 sum=vec3(0.0);\n"
        "  float weight=0.0;\n"
        "  if(center.a>0.99) {\n"
        "    sum+=center.rgb;\n"
        "    weight++;\n"
        "  }\n"
        "  if(tap1.a>0.99) {\n"
        "    sum+=tap1.rgb;\n"
        "    weight++;\n"
        "  }\n"
        "  if(tap2.a>0.99) {\n"
        "    sum+=tap2.rgb;\n"
        "    weight++;\n"
        "  }\n"
        "\n"
        " /*  sum+=center.rgb*center.a;\n"
        "     weight+=center.a;\n"
        "     sum+=tap1.rgb*tap1.a;\n"
        "     weight+=tap1.a;\n"
        "     sum+=tap2.rgb*tap2.a;\n"
        "     weight+=tap2.a;*/\n"
        "\n"
        "  if(center.a>0.9)\n"
            //we ouput the center pixel if it is not a empty edge
        "    gl_FragColor = center;\n"
        "  else\n"
             //else weighted mean of neighbourhood
        "    gl_FragColor = vec4(sum/weight,min(weight,1.0));\n"
        "}\n"
    );

    // upload the program
    m_expandSilhouetteProgram->init();
    // configure parameters
    param_expandViewTex = m_expandSilhouetteProgram->getParameterID("viewSubScatTex", false);


    /* splat rendering program */
    m_renderSplatsProgram = m_renderer->createGPUProgram("SplatSSS Splat Rendering Program");
    m_renderSplatsProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "attribute vec2 billboardOffset;\n"
        "//uniform sampler2D viewSurfPos;\n"
        "//uniform sampler2D translucencyTex;\n"
        "uniform float billboardRadius;\n"
        "\n"
        "varying vec3 splatOrigin;\n"
        "varying vec4 pixelProj;\n"
        "\n"
        "void main() {\n"
           // Vertex (splat center) position in camera space
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n"
           // Translate by current offset to generate the billboard
        "  vMV.xy += billboardOffset*billboardRadius; //*length(gl_MultiTexCoord0.xyz);\n"
           // Vertex position in clip space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n"
           // Give splat origin and screen space coordinates to fragment program
        "  splatOrigin = gl_Vertex.xyz;\n"
        "  pixelProj = vMVP;\n"
        "\n"
           //out, contains splat color
        "  gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "  gl_Position = vMVP;\n"
        "}\n"
    );

    m_renderSplatsProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "uniform sampler2D viewSurfPos;\n"
        "uniform sampler2D translucencyTex;\n"
        "uniform float billboardRadius;\n"
        "\n"
         // Corresponds to a light source in the dipole theory
        "varying vec3 splatOrigin;\n"
         // Pixel projected texture coordinate
        "varying vec4 pixelProj;\n"
        "\n"
        "void main() {\n"
           // Calculate visible surfaces
        "  vec3 visSurfPos = texture2D( viewSurfPos, 0.5 + 0.5*pixelProj.st/pixelProj.ww ).rgb;\n"
           // Distance from the splat center in screen space
        "  float dist = length(visSurfPos - splatOrigin);\n"
           // Texture coordinate in 1D diffusion profile
        "  float dist2splatCenter = dist/(billboardRadius); //*length(gl_TexCoord[0].xyz)\n"
           // gl_TexCoord[0].xyz contains the splat color at origin
        "  vec3 finalPixelColor=gl_TexCoord[0].xyz * texture2D( translucencyTex, vec2(dist2splatCenter,0.5)).rgb;\n"
        "  gl_FragColor = vec4(finalPixelColor, 0.0);\n"
        "}\n"
    );

    // upload the program
    m_renderSplatsProgram->init();
    // configure parameters
    attrib_renderSplatsBillboardOffset = m_renderSplatsProgram->getAttributeID("billboardOffset", false);
    param_renderSplatsViewSurfacePos = m_renderSplatsProgram->getParameterID("viewSurfPos", false);
    param_renderSplatsTranslucencyTex = m_renderSplatsProgram->getParameterID("translucencyTex", false);
    param_renderSplatsBillboardRadius = m_renderSplatsProgram->getParameterID("billboardRadius", false);


    /* final contribution program */
    m_finalContributionProgram = m_renderer->createGPUProgram("SplatSSS Final Rendering Program");
    m_finalContributionProgram->setSource(GPUProgram::EVertexProgram,
        "#version 120\n"
        "//uniform sampler2D subSurf;\n"
        "//uniform sampler2D albedoTex;\n"
        "uniform float sampleScale;\n"
        "uniform vec3 lightPos;\n"
        "uniform vec3 lightSpecColor;\n"
        "uniform vec3 lightDir;\n"
        "uniform float lightAperture;\n"
        "\n"
        "varying vec4 pixelProj;\n"
        "varying vec3 normal;\n"
        "varying vec3 surfPos;\n"
        "varying vec3 surfToLight;\n"
        "varying vec3 halfVec;\n"
        "\n"
        "void main() {\n"
           // Vertex position in camera space
        "  vec4 vMV = gl_ModelViewMatrix * gl_Vertex;\n"
           // Vertex position in clip space
        "  vec4 vMVP = gl_ProjectionMatrix * vMV;\n"
        "  \n"
        "  pixelProj = vMVP;\n"
        "  surfPos = gl_Vertex.xyz;\n"
        "  surfToLight = lightPos - gl_Vertex.xyz;\n"
        "  vec3 surfToView = gl_ModelViewMatrixInverse[3].xyz - gl_Vertex.xyz;\n"
           // Half vector to compute Blinn specular
        "  halfVec = ((surfToView+surfToLight)/2.0);\n"
        "  normal = gl_Normal;\n" // out
        "  gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "  gl_Position = vMVP;\n"
        "}\n"
    );

    m_finalContributionProgram->setSource(GPUProgram::EFragmentProgram,
        "#version 120\n"
        "uniform sampler2D subSurf;\n"
        "uniform sampler2D albedoTex;\n"
        "uniform float sampleScale;\n"
        "uniform vec3 lightPos;\n"
        "uniform vec3 lightSpecColor;\n"
        "uniform vec3 lightDir;\n"
        "uniform float lightAperture;\n"
        "\n"
        "varying vec4 pixelProj;\n"
        "varying vec3 normal;\n"
        "varying vec3 surfPos;\n"
        "varying vec3 surfToLight;\n"
        "varying vec3 halfVec;\n"
        "\n"
        "void main() {\n"
           // Visible surfaces
        "  \n"
          // Subsurface light
        "  vec4 subsurfaceContrib = texture2D( subSurf, (0.5 + 0.5 * pixelProj.st / pixelProj.ww) );\n"
          // Sqrt because light is musplipied two times by albedo
        "  vec3 surfaceAlbedo = sqrt(texture2D( albedoTex, gl_TexCoord[0].st ).rgb);\n"
        "  \n"
           // Compute some data for spot ligth specular"
        "  vec3 normalNorm = normalize(normal);\n"
        "  vec3 halfVecNorm = normalize(halfVec);\n"
        "  vec3 surfToLightNorm = normalize(surfToLight);\n"
        "  float specIntensity = clamp(dot(normalNorm, halfVecNorm),0.0,1.0);\n"
        "  float dot3SpotCone = dot(lightDir,-surfToLightNorm);\n"
        "  float spotCone = clamp( (dot3SpotCone-lightAperture)/(1.0-lightAperture) ,0.0,1.0);\n"
        "  float dot3lamb = dot(normal,surfToLightNorm);\n"
        "  float specContrib = clamp(dot3lamb*10.0,0.0,1.0);\n"
        "  \n"
           // Specular component
        "  vec3 specular = specIntensity * surfaceAlbedo * lightSpecColor * specContrib * spotCone * pow(specIntensity,64.0);\n"
        "  \n"
           // Final color
        "  gl_FragColor = vec4(specular + surfaceAlbedo*subsurfaceContrib.rgb/sampleScale,0.0);\n"
        "}\n"
    );

    // upload the program
    m_finalContributionProgram->init();
    // configure parameters
    param_finalContribSubSurf = m_finalContributionProgram->getParameterID("subSurf", false);
    param_finalContribAlbedoTex = m_finalContributionProgram->getParameterID("albedoTex", false);
    param_finalContribSampleScale = m_finalContributionProgram->getParameterID("sampleScale", false);
    param_finalContribLightAperture = m_finalContributionProgram->getParameterID("lightAperture", false);
    param_finalContribLightSpecColor = m_finalContributionProgram->getParameterID("lightSpecColor", false);
    param_finalContribLightDir = m_finalContributionProgram->getParameterID("lightDir", false);
    param_finalContribLightPos = m_finalContributionProgram->getParameterID("lightPos", false);


    m_initialized = true;
}

void DirectShaderManager::configure(const BSDF *bsdf,
            const Luminaire *luminaire, const SpotLight &spot,
            const Point &camPos, bool faceNormals) {
    Shader *bsdfShader = m_renderer->getShaderForResource(bsdf);
    Shader *lumShader = (luminaire == NULL) ? NULL
        : m_renderer->getShaderForResource(luminaire);
    Shader *subsurfaceShader = bsdfShader;

    std::ostringstream oss;

    if (bsdfShader == NULL || (luminaire != NULL && lumShader == NULL)) {
        /* Unsupported! */
        m_renderer->setColor(Spectrum(0.0f));
        return;
    }

   GPUProgram* program = NULL;

   //if (bsdfShader != NULL)
   //     program = m_renderer->createGPUProgram(configName);

    bool anisotropic = bsdf->getType() & BSDF::EAnisotropic;

    m_targetConfig = DirectProgramConfiguration(bsdfShader, lumShader, faceNormals);
    m_targetConfig.toString(oss);
    std::string configName = oss.str();
    std::map<std::string, ProgramAndConfiguration>::iterator it =
        m_programs.find(configName);

    if (it != m_programs.end()) {
        // A program for this configuration has been created previously
        m_current = (*it).second;
        program = m_current.program;
    } else {
        // No program for this particular combination exists -- create one
        program = m_renderer->createGPUProgram(configName);

		if (faceNormals) {
			/* Generate face normals in a geometry shader */
	
    		if (!m_renderer->getCapabilities()->isSupported(
					RendererCapabilities::EGeometryShaders))
				Log(EError, "Face normals require geometry shader support!");
			if (anisotropic)
				Log(EError, "Anisotropy and face normals can't be combined at the moment");
	
			oss.str("");
			oss << "#version 120" << endl
				<< "#extension GL_EXT_geometry_shader4 : enable" << endl
				<< "varying in vec3 lightVec_vertex[3], camVec_vertex[3];" << endl
				<< "varying in vec2 uv_vertex[3];" << endl
				<< "varying in vec3 vertexColor_vertex[3];" << endl
				<< "varying out vec3 normal;" << endl
				<< "varying out vec3 lightVec, camVec;" << endl
				<< "varying out vec2 uv;" << endl
				<< "varying out vec3 vertexColor;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   vec3 p0 = gl_PositionIn[0].xyz / gl_PositionIn[0].w;" << endl
				<< "   vec3 p1 = gl_PositionIn[1].xyz / gl_PositionIn[1].w;" << endl
				<< "   vec3 p2 = gl_PositionIn[2].xyz / gl_PositionIn[2].w;" << endl
				<< "   normal = normalize(cross(p1 - p0, p2 - p0));" << endl
				<< "   gl_Position = vec4(0.0);" << endl
				<< "   lightVec = camVec = vec3(0.0);" << endl
				<< "   for (int i=0; i<gl_VerticesIn; ++i) {" << endl
				<< "      gl_Position = gl_PositionIn[i];" << endl
				<< "      uv = uv_vertex[i];" << endl
				<< "      vertexColor = vertexColor_vertex[i];" << endl
				<< "      lightVec = lightVec_vertex[i];" << endl
				<< "      camVec = camVec_vertex[i];" << endl
				<< "      EmitVertex();" << endl
				<< "   }" << endl
				<< "   EndPrimitive();" << endl
				<< "}" << endl;

			program->setMaxVertices(3); 
			program->setSource(GPUProgram::EGeometryProgram, oss.str());
		}

		/* Vertex program */
		oss.str("");
        oss << "#version 120" << endl;
		if (anisotropic)
			oss << "varying vec3 tangent;" << endl;
		oss << "uniform vec3 lightPos, camPos;" << endl;

		if (!faceNormals) {
			oss << "varying vec3 lightVec, camVec;" << endl
				<< "varying vec2 uv;" << endl
				<< "varying vec3 normal;" << endl
				<< "varying vec3 vertexColor;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   uv = gl_MultiTexCoord0.xy;" << endl
				<< "   camVec = camPos - gl_Vertex.xyz;" << endl
				<< "   lightVec = lightPos - gl_Vertex.xyz;" << endl
				<< "   gl_Position = ftransform();" << endl
				<< "   vertexColor = gl_Color.rgb;" << endl
				<< "   normal = gl_Normal;" << endl;
		} else {
			oss << "varying vec3 lightVec_vertex, camVec_vertex;" << endl
				<< "varying vec2 uv_vertex;" << endl
				<< "varying vec3 vertexColor_vertex;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   uv_vertex = gl_MultiTexCoord0.xy;" << endl
				<< "   camVec_vertex = camPos - gl_Vertex.xyz;" << endl
				<< "   lightVec_vertex = lightPos - gl_Vertex.xyz;" << endl
				<< "   gl_Position = ftransform();" << endl
				<< "   vertexColor_vertex = gl_Color.rgb;" << endl;
		}
		if (anisotropic)
			oss << "   tangent = gl_MultiTexCoord1.xyz;" << endl;
		oss << "}" << endl;

		program->setSource(GPUProgram::EVertexProgram, oss.str());
		oss.str("");

		oss << "#version 120" << endl
			<< endl
			<< "/* Uniform inputs */" << endl
			<< "uniform samplerCube shadowMap;" << endl
            << "uniform vec3 lightPower;" << endl
			<< "uniform float nearClip, invClipRange, minDist;" << endl
			<< "uniform bool diffuseSources, diffuseReceivers;" << endl
			<< "varying vec3 vertexColor;" << endl
			<< endl
			<< "/* Inputs <- Vertex program */" << endl
			<< "varying vec3 normal, lightVec, camVec;" << endl
			<< "varying vec2 uv;" << endl;
		if (anisotropic)
			oss << "varying vec3 tangent;" << endl;

		oss << endl
			<< "/* Some helper functions for BSDF implementations */" << endl
			<< "float cosTheta(vec3 v) { return v.z; }" << endl
			<< "float sinTheta2(vec3 v) { return 1.0-v.z*v.z; }" << endl
			<< "float sinTheta(vec3 v) { float st2 = sinTheta2(v); if (st2 <= 0) return 0.0; else return sqrt(sinTheta2(v)); }" << endl
			<< "float tanTheta(vec3 v) { return sinTheta(v)/cosTheta(v); }" << endl
			<< endl;

		std::string bsdfEvalName, lumEvalName;
		m_targetConfig.generateCode(oss, bsdfEvalName, lumEvalName);

		oss << "void main() {" << endl
			<< "   /* Set up an ONB */" << endl
			<< "   vec3 N = normalize(normal);" << endl;
		if (anisotropic) {
			oss << "   vec3 S = normalize(tangent - dot(tangent, N)*N);" << endl;
		} else {
			oss << "   vec3 S;" << endl
				<< "   if (abs(N.x) > abs(N.y)) {" << endl
				<< "      float invLen = 1.0 / sqrt(N.x*N.x + N.z*N.z);" << endl
				<< "      S = vec3(-N.z * invLen, 0.0, N.x * invLen);" << endl
				<< "   } else {" << endl
				<< "      float invLen = 1.0 / sqrt(N.y*N.y + N.z*N.z);" << endl
				<< "      S = vec3(0.0, -N.z * invLen, N.y * invLen);" << endl
				<< "   }" << endl;
		}
		oss << "   vec3 T = cross(N, S);" << endl
			<< endl
			<< "   /* Compute shadows */" << endl
			<< "   float d = length(lightVec);" << endl
			<< "   vec3 nLightVec = lightVec/d, absLightVec = abs(lightVec);" << endl
			<< "   float depth = max(max(absLightVec.x, absLightVec.y), absLightVec.z);" << endl
			<< "   depth = (depth-nearClip) * invClipRange - 0.005;" << endl
			<< "   float shadow = textureCube(shadowMap, nLightVec).r > depth ? 1.0 : 0.0;" << endl
			<< endl
			<< "   /* Shading */" << endl
			<< "   vec3 nCamVec = normalize(camVec);" << endl
			<< "   vec3 wo = vec3(dot(S, nLightVec)," << endl
			<< "                  dot(T, nLightVec)," << endl
			<< "                  dot(N, nLightVec));" << endl
			<< "   vec3 wi = vec3(dot(S, nCamVec)," << endl
			<< "                  dot(T, nCamVec)," << endl
			<< "                  dot(N, nCamVec));" << endl
			<< "   vec3 contrib = lightPower;" << endl;
			//<< "   if (!diffuseSources)" << endl 
			//<< "      contrib *= " << vplEvalName;
			//if (vpl.type == ESurfaceVPL)
			//	oss << "(vplUV, vplWi, vplWo);" << endl;
			//else
			//	oss << "_dir(vplWo);" << endl;
		oss << "   if (d < minDist) d = minDist;" << endl
			<< "   if (!diffuseReceivers)" << endl
			<< "      contrib *= "<< bsdfEvalName << "(uv, wi, wo);" << endl
			<< "   else" << endl
			<< "      contrib *= " << bsdfEvalName << "_diffuse(uv, wi, wo);" << endl
			<< "   gl_FragColor.rgb = contrib";
		//if (vpl.type == ESurfaceVPL || (vpl.type == ELuminaireVPL 
		//		&& (vpl.luminaire->getType() & Luminaire::EOnSurface)))
		//	oss << " * (shadow * abs(cosTheta(wo) * cosTheta(vplWo)) / (d*d))";
		//else 
		oss << " * (shadow * abs(cosTheta(wo)) / (d*d))";
		if (luminaire != NULL) {
			oss << endl;
			oss << "                      + " << lumEvalName << "_area(uv)"
				<< " * " << lumEvalName << "_dir(wi);" << endl;
		} else {
			oss << ";" << endl;
		}
		oss << "   gl_FragColor.a = 1.0;" << endl
			<< "}" << endl;

		program->setSource(GPUProgram::EFragmentProgram, oss.str());
		try {
			program->init();
		} catch (const std::exception &) {
			Log(EWarn, "Unable to compile the following VPL program:\n%s", oss.str().c_str());
			throw;
		}

		m_targetConfig.resolve(program);
		m_targetConfig.param_shadowMap = program->getParameterID("shadowMap", false);
		m_targetConfig.param_camPos = program->getParameterID("camPos", false);
		m_targetConfig.param_lightPower = program->getParameterID("lightPower", false);
		m_targetConfig.param_lightPos = program->getParameterID("lightPos", false);
		m_targetConfig.param_nearClip = program->getParameterID("nearClip", false);
		m_targetConfig.param_invClipRange = program->getParameterID("invClipRange", false);
		m_targetConfig.param_minDist = program->getParameterID("minDist", false);
		m_targetConfig.param_diffuseSources = program->getParameterID("diffuseSources", false);
		m_targetConfig.param_diffuseReceivers = program->getParameterID("diffuseReceivers", false);
		m_current.program = program;
		m_current.config = m_targetConfig;
		m_programs[configName] = m_current;
		program->incRef();

        // ToDo: Shadow program
    }

    program->bind();
    //m_shadowMap->bind();

	const DirectProgramConfiguration &config = m_current.config;

	//program->setParameter(config.param_shadowMap, m_shadowMap);
	program->setParameter(config.param_lightPos, spot.pos);
	program->setParameter(config.param_camPos, camPos);
	//program->setParameter(config.param_diffuseSources, m_diffuseSources);
	//Spectrum power = spot.color;
	//if (m_diffuseSources && vpl.type == ESurfaceVPL)
	//	power *= vpl.its.shape->getBSDF()->getDiffuseReflectance(vpl.its) * INV_PI;
    // ToDo: 100 is just a migic number, the light is too dark without it
	program->setParameter(config.param_lightPower, spot.color * 100 /* power */);
	program->setParameter(config.param_diffuseReceivers, m_diffuseReceivers);
	program->setParameter(config.param_nearClip, m_nearClip);
	program->setParameter(config.param_invClipRange, m_invClipRange);
	program->setParameter(config.param_minDist, m_minDist);

	int textureUnitOffset = 1;
	m_targetConfig.bind(program, config, textureUnitOffset);
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
    if (m_current.program && m_current.program->isBound()) {
        m_targetConfig.unbind();
        m_current.program->unbind();
        //m_shadowMap->unbind();
    }
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

