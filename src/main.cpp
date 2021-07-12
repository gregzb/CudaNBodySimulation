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
void process_input(GLFWwindow *window);

GLFWwindow* initialize_glfw() {
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

void init_sphere_buffers(unsigned int &VBO, unsigned int &VAO, unsigned int &EBO, const std::vector<vertex> &vertices, const std::vector<unsigned int> &indices) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*vertices.size(), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*indices.size(), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)(sizeof(float)*3));
    glEnableVertexAttribArray(1);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0); 

    // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    glBindVertexArray(0);
}

void enable_gl_settings() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_MULTISAMPLE);

    // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();

    GLFWwindow* window;
    if (!(window=initialize_glfw())) {
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

    shader<GL_VERTEX_SHADER> test_vertex_shader("shaders/body_vertex.glsl");
    shader<GL_FRAGMENT_SHADER> test_fragment_shader("shaders/body_fragment.glsl");
    shader_program test_shader({test_vertex_shader.get_id(), test_fragment_shader.get_id()});
    test_vertex_shader.release();
    test_fragment_shader.release();

    test_shader.add_uniform("model");
    test_shader.add_uniform("normal_model");
    test_shader.add_uniform("view_projection");
    test_shader.add_uniform("light_pos");
    test_shader.add_uniform("view_pos");

    std::vector<vertex> vertices;
    std::vector<unsigned int> indices;

    generateSphereMesh(vertices, indices, 16, 32);

    unsigned int VBO, VAO, EBO;
    init_sphere_buffers(VBO, VAO, EBO, vertices, indices);

    unsigned int buffer;
    glGenBuffers(1, &buffer);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    std::vector<float> dat{2, 0, -3, 1, 0, 0, -2, 0, -3, 0, 0, 1};
    glBufferData(GL_ARRAY_BUFFER, dat.size() * sizeof(float), dat.data(), GL_STATIC_DRAW);

    glBindVertexArray(VAO);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)0);
    glVertexAttribDivisor(2, 1);

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)(sizeof(float)*3));
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(0);

    enable_gl_settings();

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
        process_input(window);

        // render
        // ------
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        test_shader.use();

        glm::mat4 model(1.0f);
        // model = glm::scale(model, {0.1f, 0.1f, 0.1f});
        // model = glm::rotate(model, (float)newTime, {0, 1, 0});
        // cam.set_target({3*std::sin(newTime), 0, -3});
        // cam.set_pos({5*std::sin(newTime), 0, 5*std::cos(newTime)});
        glm::mat4 vp_matrix = cam.get_projection_matrix() * cam.get_view_matrix();

        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         std::cout << mvp[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // break;

        glUniformMatrix4fv(test_shader.get_uniform_location("view_projection"), 1, GL_FALSE, &vp_matrix[0][0]);
        glUniformMatrix4fv(test_shader.get_uniform_location("model"), 1, GL_FALSE, &model[0][0]);
        glm::mat3 normal_model_view = glm::transpose(glm::inverse(glm::mat3(model)));
        glUniformMatrix3fv(test_shader.get_uniform_location("normal_model"), 1, GL_FALSE, &normal_model_view[0][0]);

        glm::vec3 light_pos {(float)std::sin(newTime)*5, 0, 2};
        glUniform3f(test_shader.get_uniform_location("light_pos"), light_pos.x, light_pos.y, light_pos.z);

        glUniform3f(test_shader.get_uniform_location("view_pos"), cam.get_pos().x, cam.get_pos().y, cam.get_pos().z);

        glBindVertexArray(VAO);
        glDrawElementsInstanced(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (void*) 0, 2);

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
void process_input(GLFWwindow *window)
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