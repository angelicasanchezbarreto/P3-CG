#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//#include <learnopengl/filesystem.h>
#include "shader_m.h"
#include "camera.h"
#include "Objeto.h"
#include "glut_ply.h"
#include <iostream>
#include "DTriangulation.h"
using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 40.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;	
float lastFrame = 0.0f;
float tiempoInicial = 0.0f, tiempoTranscurrido = 0.0f;

// lighting
glm::vec3 lightPos(1.2f, 30.0f, 2.0f);

Esfera esfera(vec3(0),2., 20, 20);
Esfera *pEsfera = new Esfera(vec3(0),2, 50, 50);
Model_PLY modelo;
vector<Objeto*> objetos;
bool boton_presionado = false;

int main() {
    //char *archivo = "../models/bunny.ply";
    //modelo.Load(archivo);

    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)     {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))     {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader program
    Shader lightingShader("../2.2.basic_lighting.vs", "../2.2.basic_lighting.fs");
    Shader lightCubeShader("../2.2.light_cube.vs", "../2.2.light_cube.fs");

    vector<vec3> puntos;
    float radius = 2.0f;
    float height = 10.0f;
    int numSegments = 20;
    int numStacks = 10;
    float angleIncrement = 2 * M_PI / numSegments;
    float heightStep = height / numStacks;

    // Generate the points for a tunnel
    for (int stack = 0; stack <= numStacks; ++stack) {
        float currentHeight = stack * heightStep;
        for (int segment = 0; segment <= numSegments; ++segment) {
            float currentAngle = segment * angleIncrement;
            float x = radius * cos(currentAngle);
            float y = radius * sin(currentAngle);
            float z = currentHeight;

            puntos.push_back({x, y, z});
        }
    }
    DTriangulation dt;
    vector<Triangle3D> triangles = dt.performTriangulation(puntos);

    GLuint vao_puntos;
    GLint POSITION_ATTRIBUTE=0;
    glGenVertexArrays( 1, &vao_puntos );
    glBindVertexArray( vao_puntos );
    GLuint vbos_puntos[1];
    glGenBuffers( 1, vbos_puntos );
    glBindBuffer( GL_ARRAY_BUFFER, vbos_puntos[0] );
    glBufferData( GL_ARRAY_BUFFER, puntos.size() * sizeof(vec3), puntos.data(), GL_STATIC_DRAW );

    /*GLuint vao_triangles;
    glGenVertexArrays(1, &vao_triangles);
    glBindVertexArray(vao_triangles);
    GLuint vbo_triangles[0];
    glGenBuffers(1, vbo_triangles);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_triangles[0]);
    glBufferData(GL_ARRAY_BUFFER, triangles.size() * sizeof(Triangle3D), triangles.data(), GL_STATIC_DRAW);*/

    // Set up the vertex attribute pointers
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), nullptr);

    glEnableVertexAttribArray( POSITION_ATTRIBUTE );
    glVertexAttribPointer( POSITION_ATTRIBUTE, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glBindVertexArray( 0 );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );

    // render loop
    while (!glfwWindowShouldClose(window)) {
        // per-frame time logic
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        // input
        processInput(window);

        // render
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShader.use();
        lightingShader.setVec3("objectColor", 1.0f, 0.5f, 0.31f);
        lightingShader.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
        lightingShader.setVec3("lightPos", lightPos);
        lightingShader.setVec3("viewPos", camera.Position);

        // view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        lightingShader.setMat4("projection", projection);
        lightingShader.setMat4("view", view);

        // world transformation
        glm::mat4 model = glm::mat4(1.0f);
        lightingShader.setMat4("model", model);

        lightCubeShader.use();
        lightCubeShader.setMat4("projection", projection);
        lightCubeShader.setMat4("view", view);
        lightCubeShader.setMat4("model", model);

        /*
        // Drawing the triangles
        for (const auto& triangle : triangles) {
            glVertex3f(triangle.p1.x, triangle.p1.y, triangle.p1.z);
            glVertex3f(triangle.p2.x, triangle.p2.y, triangle.p2.z);
            glVertex3f(triangle.p3.x, triangle.p3.y, triangle.p3.z);
        }*/

        glBindVertexArray(vao_puntos);
        glDrawArrays(GL_POINTS,0, puntos.size());
        glBindVertexArray(0);

        /*glBindVertexArray(vao_triangles);
        glDrawArrays(GL_TRIANGLES, 0, triangles.size() * 3);
        glBindVertexArray(0);*/

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
 //   glDeleteVertexArrays(1, &cubeVAO);
 //   glDeleteVertexArrays(1, &lightCubeVAO);
 //   glDeleteBuffers(1, &VBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
void processInput(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS){
        boton_presionado = true;
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_RELEASE){
        if (boton_presionado) {
            float x = rand()%10;
            float y = rand()%10;
            float z = rand()%10;
            Objeto *pE = new Esfera(glm::vec3(x,y,z));
            pE->centro = vec3(x,y,z);
            pE->v0 = 20;
            pE->a0 = 50 + rand() % 20;
            pE->x0 = x;
            pE->y0 = y;
            pE->vao = esfera.vao;
            pE->indices_size = esfera.indices_size;
            pE->tiempo_inicial = static_cast<float>(glfwGetTime());
            objetos.emplace_back(pE);
            boton_presionado = false;
            //cout << endl << x << " " << y << " " << z << " " << pE->a0;

        }
    }

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;
    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}
