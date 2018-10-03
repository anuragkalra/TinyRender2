/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include "core/renderpass.h"
#include <bsdfs/diffuse.h>
#include <string>
#include <fstream>

TR_NAMESPACE_BEGIN

bool RenderPass::init(const Config& config) {
    fs::path file("../../../../src/shaders/");
    file = (config.tomlFile.parent_path() / file).make_preferred();
    shadersFilePath = file.string();

    // Create matrices
    modelMat = glm::scale(glm::mat4(1.f), glm::vec3(1.f));
    normalMat = glm::transpose(glm::inverse(modelMat));

    // Camera
    camera.SetPosition(config.camera.o);
    camera.SetLookAt(config.camera.at);
    camera.SetUp(config.camera.up);
    camera.SetClipping(.01, 1000);
    camera.SetFOV(glm::radians(config.camera.fov));
    camera.SetViewport(0, 0, config.width, config.height);
    camera.camera_scale = .005f;
    camera.max_pitch_rate = 0.005f;
    camera.max_heading_rate = 0.005f;

    // Light source
    size_t idxEmitter = getFirstLight();
    if (idxEmitter != -1) {
        lightPos = scene.worldData.shapesCenter[idxEmitter];
        lightPower = scene.emitters[0].getPower() / (2 * M_PI); // Convert area light to point light
    } else {
        std::cout << "Can't find light source!" << std::endl;
    }

    return true;
}

void RenderPass::cleanUp() {
    SDL_GL_DeleteContext(contextGL);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void RenderPass::render() {

    if (!isSaved) { // save first frame
        std::unique_ptr<GLfloat> data(new GLfloat[3 * nPixel]);
        glReadPixels(0, 0, width, height, GL_RGB, GL_FLOAT, data.get());
        isSaved = save(data.get());
    }
}

bool RenderPass::save(GLfloat* data) {
    fs::path p = scene.config.tomlFile;
    fs::path newP = p.parent_path() / fs::path(p.stem().string() + p.extension().string());
    saveEXR(data, newP.replace_extension("exr").string(), scene.config.width, scene.config.height);
    return true;
}

bool RenderPass::initOpenGL(int width, int height) {
    this->width = width;
    this->height = height;
    nPixel = width * height;

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cout << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return false;
    }

    // Create window
    window = SDL_CreateWindow("TinyRender (Real-Time)", 100, 100, width, height, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if (!window) {
        std::cout << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return false;
    }

    // Use OpenGL core profile (disable legacy stuff, it ain't the 90's anymore...)
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetSwapInterval(1); // enable vsync

    // Create OpenGL context
    contextGL = SDL_GL_CreateContext(window);
    if (contextGL == nullptr) {
        std::cout << "SDL_GL_CreateContext error: " << SDL_GetError() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Init GLEW (needs to be called just after creating the OpenGL context)
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        std::cout << "glewInit error: " << glewGetErrorString(err) << std::endl;
    }

    //clear and swap
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    SDL_GL_SwapWindow(window);

    return true;
}

GLuint RenderPass::compileProgram(GLuint vs, GLuint fs) {
    printf("Linking program...\n");

    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    // Check for errors
    GLint isLinked = GL_FALSE;
    glGetProgramiv(program, GL_LINK_STATUS, &isLinked);
    if (isLinked == GL_FALSE) {
        GLint maxLength = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

        std::vector<GLchar> errorLog(maxLength);
        glGetProgramInfoLog(program, maxLength, &maxLength, &errorLog[0]);

        glDeleteProgram(program);
        glDeleteShader(vs);
        glDeleteShader(fs);

        printf("%s\n", &errorLog[0]);
        exit(EXIT_FAILURE);
    }

    glDetachShader(program, vs);
    glDetachShader(program, fs);

    return program;
}

GLuint RenderPass::compileShader(const char* shaderPath_, GLenum shaderType) {
    std::string shaderPath = shadersFilePath;
    shaderPath += shaderPath_;

    printf("Compiling shader: %s\n", shaderPath.c_str());

    std::string code = readFile(shaderPath.c_str());
    GLuint id = glCreateShader(shaderType);
    char const* codePtr = code.c_str();
    glShaderSource(id, 1, &codePtr, NULL);
    glCompileShader(id);

    // Check for errors
    GLint isCompiled = GL_FALSE;
    glGetShaderiv(id, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        int maxLength = 0;
        glGetShaderiv(id, GL_INFO_LOG_LENGTH, &maxLength);

        std::vector<GLchar> errorLog(maxLength);
        glGetShaderInfoLog(id, maxLength, &maxLength, &errorLog[0]);

        glDeleteShader(id);

        printf("%s\n", &errorLog[0]);
        exit(EXIT_FAILURE);
    }

    return id;
}

std::string RenderPass::readFile(const char* filePath) {
    std::string str;

    std::ifstream file(filePath, std::ios::in);
    if (file.is_open()) {
        std::stringstream sstr;
        sstr << file.rdbuf();
        str = sstr.str();
        file.close();
    } else {
        printf("Error: Can't read file '%s'\n", filePath);
        exit(EXIT_FAILURE);
    }

    return str;
}

void RenderPass::buildVBO(size_t objectIdx) {
    const tinyobj::attrib_t& sa = scene.worldData.attrib;
    const tinyobj::shape_t& s = scene.worldData.shapes[objectIdx];

    GLObject& obj = objects[objectIdx];

    obj.nVerts = s.mesh.indices.size();
    obj.vertices.resize(obj.nVerts * N_ATTR_PER_VERT);
    int k = 0;
    for (size_t i = 0; i < s.mesh.indices.size(); i++) {
        int idx = s.mesh.indices[i].vertex_index;
        int idx_n = s.mesh.indices[i].normal_index;

        // Position
        obj.vertices[k + 0] = sa.vertices[3 * idx + 0];
        obj.vertices[k + 1] = sa.vertices[3 * idx + 1];
        obj.vertices[k + 2] = sa.vertices[3 * idx + 2];

        // Normal
        float nx = sa.normals[3 * idx_n + 0];
        float ny = sa.normals[3 * idx_n + 1];
        float nz = sa.normals[3 * idx_n + 2];
        float norm = std::sqrt(nx * nx + ny * ny + nz * nz);
        obj.vertices[k + 3] = nx / norm;
        obj.vertices[k + 4] = ny / norm;
        obj.vertices[k + 5] = nz / norm;

        k += N_ATTR_PER_VERT;
    }

    // VBO
    glGenVertexArrays(1, &obj.vao);
    glBindVertexArray(obj.vao);

    glGenBuffers(1, &obj.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
                 (GLvoid*) (&obj.vertices[0]),
                 GL_STATIC_DRAW);
}

void RenderPass::buildVAO(size_t objectIdx) {
    glBindVertexArray(objects[objectIdx].vao);
    glBindBuffer(GL_ARRAY_BUFFER, objects[objectIdx].vbo);

    // Define vertices attributes
    glEnableVertexAttribArray(posAttrib);
    glEnableVertexAttribArray(normalAttrib);
    glVertexAttribPointer(posAttrib,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          sizeof(GLfloat) * N_ATTR_PER_VERT,
                          (GLvoid*) (0 * sizeof(GLfloat)));
    glVertexAttribPointer(normalAttrib,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          sizeof(GLfloat) * N_ATTR_PER_VERT,
                          (GLvoid*) (3 * sizeof(GLfloat)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void RenderPass::assignShader(GLObject& obj,
                              const tinyobj::shape_t& s,
                              const std::vector<std::unique_ptr<BSDF>>& bsdfs) {
}

size_t RenderPass::getFirstLight() {
    const auto& shapes = scene.worldData.shapes;
    for (size_t i = 0; i < shapes.size(); i++) {
        const tinyobj::shape_t& s = shapes[i];
        int materialId = s.mesh.material_ids[0];
        const BSDF* bsdf = scene.bsdfs[materialId].get();
        if (bsdf->isEmissive()) {
            return i;
        }
    }
    return -1;
}

void RenderPass::updateCamera(SDL_Event& e) {
    switch (e.type) {

        case SDL_KEYDOWN:
            if (e.key.keysym.sym == SDLK_w) {
                camera.Move(FORWARD);
            } else if (e.key.keysym.sym == SDLK_a) {
                camera.Move(LEFT);
            } else if (e.key.keysym.sym == SDLK_s) {
                camera.Move(BACK);
            } else if (e.key.keysym.sym == SDLK_d) {
                camera.Move(RIGHT);
            }
            break;

        case SDL_MOUSEBUTTONDOWN:camera.move_camera = true;
            break;

        case SDL_MOUSEBUTTONUP:camera.move_camera = false;
            break;

        case SDL_MOUSEMOTION: {
            /* If the mouse is moving to the left */
            if (e.motion.xrel < 0)
                camera.Move2D(e.motion.x, e.motion.y);
            else if (e.motion.xrel > 0)
                camera.Move2D(e.motion.x, e.motion.y);
        }
            break;
    }
}

TR_NAMESPACE_END
