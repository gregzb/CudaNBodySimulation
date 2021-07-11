#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include "shader.hpp"
#include "shader_program.hpp"

#include "shapes.hpp"
#include "vertex.hpp"
#include "utils.hpp"
#include "camera.hpp"
#include "graphics_settings.hpp"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

GLFWwindow* initializeGLFW() {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "N-Body Simulation", NULL, NULL);
    if (window == NULL)
    {
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    return window;
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();

    GLFWwindow* window;
    if (!(window=initializeGLFW())) {
        std::cout << "Failed to create GLFW window" << std::endl;
        return -1;
    }

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // shader<GL_VERTEX_SHADER> circle_vertex_shader("shaders/circle_vertex.glsl");
    // shader<GL_FRAGMENT_SHADER> circle_fragment_shader("shaders/circle_fragment.glsl");
    // shader_program circle_shader({circle_vertex_shader.get_id(), circle_fragment_shader.get_id()});
    // circle_vertex_shader.release();
    // circle_fragment_shader.release();

    shader<GL_VERTEX_SHADER> test_vertex_shader("shaders/test_shader_vertex.glsl");
    shader<GL_FRAGMENT_SHADER> test_fragment_shader("shaders/test_shader_fragment.glsl");
    shader_program test_shader({test_vertex_shader.get_id(), test_fragment_shader.get_id()});
    test_vertex_shader.release();
    test_fragment_shader.release();

    test_shader.add_uniform("mvp");

    std::vector<vertex> vertices;
    std::vector<unsigned int> indices;

    generateSphereMesh(vertices, indices, 10, 10);

    // std::vector<float> vertices {
    //     0.5f,  0.5f, 0.0f,  // top right
    //     0.5f, -0.5f, 0.0f,  // bottom right
    //     -0.5f, -0.5f, 0.0f,  // bottom left
    //     -0.5f,  0.5f, 0.0f   // top left 
    // };

    // std::vector<vertex> vertices {
    //     {{0.5f, 0.5f, 0.0f}},
    //     {{0.5f, -0.5f, 0.0f}},
    //     {{-0.5f, -0.5f, 0.0f}},
    //     {{-0.5f, 0.5f, 0.0f}}
    // };

    // std::vector<unsigned int> indices {  // note that we start from 0!
    //     0, 3, 1,   // first triangle
    //     1, 3, 2    // second triangle
    // };

    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*vertices.size(), vertices.data(), GL_STATIC_DRAW);

    // for (unsigned i = 0; i < vertices.size(); i++) {
    //     std::cout << vertices[i].get_pos() << " ";
    // }
    // std::cout << std::endl;

    // for (unsigned i = 0; i < indices.size(); i++) {
    //     std::cout << indices[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << sizeof(vertex) << std::endl;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*indices.size(), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0); 

    // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_MULTISAMPLE);

    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    // glm::mat4 view          = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    // glm::mat4 projection    = glm::perspective(glm::radians(75.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.01f, 1000.0f);
    // view       = glm::translate(view, glm::vec3(0.0f, 0.0f, -3.0f));

    camera cam(glm::radians(75.0f), {0, 0, 0}, {0, 0, -1});

    double lastTime = glfwGetTime();

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        double newTime = glfwGetTime();
        double dt = newTime-lastTime;
        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        test_shader.use();

        glm::mat4 model(1.0f);
        model = glm::translate(model, {3*std::sin(newTime), 0, -3});
        cam.set_target({3*std::sin(newTime), 0, -3});
        cam.set_pos({5*std::sin(newTime), 0, 5*std::cos(newTime)});
        glm::mat4 vp_matrix = cam.get_vp_matrix();
        glm::mat4 mvp = vp_matrix * model;

        glUniformMatrix4fv(test_shader.get_uniform_location("mvp"), 1, GL_FALSE, &mvp[0][0]);

        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (void*) 0);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();

        lastTime = newTime;
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}