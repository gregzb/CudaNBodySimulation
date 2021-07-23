#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <random>
#include <cmath>

#include <glm/gtc/random.hpp>
#include <glm/gtx/norm.hpp>

#include "shader.hpp"
#include "shader_program.hpp"

#include "shapes.hpp"
#include "vertex.hpp"
#include "utils.hpp"
#include "camera.hpp"
#include "graphics_settings.hpp"
#include "body.hpp"
#include "nbody_simulation.hpp"

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

GLFWwindow* initialize_glfw() {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_FALSE); // single buffering, inherently disables vsync?

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

void initBodies(std::vector<body> &bodies) {
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen(0); // keep it reproducible for now
    std::uniform_real_distribution<float> radius_distrib(1, 2);
    std::uniform_real_distribution<float> phi_distrib(0, (float)M_PI);
    std::uniform_real_distribution<float> theta_distrib(0, (float)M_PI*2);

    glm::vec3 center(0, 0, 0);
    bodies.push_back({center, {0, 0, 0}, 10000000});

    float initial_velocity_mag = 0.02;

    for (int i = 0; i < 100; i++) {
        float r = radius_distrib(gen);
        float phi = phi_distrib(gen);
        float theta = theta_distrib(gen);
        glm::vec3 pos(r*std::cos(theta)*std::sin(phi), r*std::sin(theta)*std::sin(phi), r*std::cos(phi));
        glm::vec3 normal(glm::normalize(pos-center));

        float epsilon = 0.01;

        glm::vec3 another;
        while(glm::length2(cross(another=glm::sphericalRand(1.0f), normal)) < epsilon);

        glm::vec3 vector_on_plane = glm::normalize(cross(normal, another)) * initial_velocity_mag;

        bodies.push_back({pos, vector_on_plane, 0.01});
    }
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();

    GLFWwindow* window;
    if (!(window = initialize_glfw())) {
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

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");

    shader<GL_VERTEX_SHADER> body_vertex_shader("shaders/body_vertex.glsl");
    shader<GL_FRAGMENT_SHADER> body_fragment_shader("shaders/body_fragment.glsl");
    shader_program body_shader({body_vertex_shader.get_id(), body_fragment_shader.get_id()});
    body_vertex_shader.release();
    body_fragment_shader.release();

    body_shader.add_uniform("model");
    body_shader.add_uniform("view");
    body_shader.add_uniform("projection");
    body_shader.add_uniform("normal_model_view");
    body_shader.add_uniform("light_pos");

    std::vector<vertex> vertices;
    std::vector<unsigned int> indices;

    generateSphereMesh(vertices, indices, 8, 8);

    unsigned int VBO, VAO, EBO;
    init_sphere_buffers(VBO, VAO, EBO, vertices, indices);

    std::vector<body> bodies;
    initBodies(bodies);

    nbody_simulation simulation(bodies, 0.1f);
    simulation.init();

    unsigned int instance_buffer;
    glGenBuffers(1, &instance_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, instance_buffer);
    std::vector<float> dat;
    for(const auto &body : bodies) {
        dat.push_back(body.pos.x);
        dat.push_back(body.pos.y);
        dat.push_back(body.pos.z);
        dat.push_back(0.2);
        dat.push_back(0.4);
        dat.push_back(1);
    }
    glBufferData(GL_ARRAY_BUFFER, dat.size() * sizeof(float), dat.data(), GL_DYNAMIC_DRAW);

    glBindVertexArray(VAO);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)0);
    glVertexAttribDivisor(2, 1);

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)(sizeof(float)*3));
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(0);

    enable_gl_settings();

    camera cam(glm::radians(60.0f), {0, 0, 0}, {0, 0, 0});

    double lastTime = glfwGetTime();

    double frame_time = 0;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        double newTime = glfwGetTime();
        double dt = newTime-lastTime;
        while(dt < (1/900.0)) {
            newTime = glfwGetTime();
            dt = newTime-lastTime;
        }
        // std::cout << dt << std::endl;
        // input
        // -----

        glfwPollEvents();
        // glfwWaitEventsTimeout(0.015);
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
        {
            static float time_step = 1.0f;
            static ImVec4 highlight_color(0.2f, 0.9f, 0.2f, 1.0f);

            ImGui::Begin("Control Panel");                          // Create a window called "Hello, world!" and append into it.

            ImGui::Text("Use this panel to control properties of the simulation");               // Display some text (you can use a format strings too)
            // ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
            // ImGui::Checkbox("Another Window", &show_another_window);

            ImGui::SliderFloat("Time Step", &time_step, 0.1f, 10.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
            ImGui::ColorEdit3("Highlight Color", (float*)&highlight_color); // Edit 3 floats representing a color

            if (ImGui::Button("Highlight Random Body")) {
                std::cout << "highlighting!" << std::endl;
            }
            //     counter++;
            // ImGui::SameLine();
            // ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        process_input(window);

        // render
        // ------
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        body_shader.use();

        // cam.set_pos({5.5*std::sin(newTime), 0, 5.5*std::cos(newTime)});
        cam.set_pos({0, 0, 5.5});

        glm::mat4 model(1.0f);
        model = glm::scale(model, {0.03f, 0.03f, 0.03f});
        const glm::mat4 &view = cam.get_view_matrix();
        const glm::mat4 &projection = cam.get_projection_matrix();

        body_shader.set_mvp(model, view, projection);
        glm::mat3 normal_model_view = glm::transpose(glm::inverse(glm::mat3(cam.get_view_matrix() * model)));
        glUniformMatrix3fv(body_shader.get_uniform_location("normal_model_view"), 1, GL_FALSE, &normal_model_view[0][0]);

        glm::vec3 light_world_pos {0, 8, 8};
        glm::vec3 light_pos = glm::vec3(view * glm::vec4(light_world_pos, 1.0f));
        glUniform3f(body_shader.get_uniform_location("light_pos"), light_pos.x, light_pos.y, light_pos.z);

        glBindVertexArray(VAO);

        simulation.step();

        std::cout << simulation.get_energy() << std::endl;

        glm::vec3 red(1, 0.2, 0.2);
        glm::vec3 blue(0.2, 0.4, 1);

        glBindBuffer(GL_ARRAY_BUFFER, instance_buffer);
        std::vector<float> dat;
        for(const auto &body : simulation.get_bodies()) {
            dat.push_back(body.pos.x);
            dat.push_back(body.pos.y);
            dat.push_back(body.pos.z);
            float mag = glm::length(body.vel);
            float mixer = glm::clamp(mag*25, 0.0f, 1.0f);
            glm::vec3 color(glm::mix(blue, red, mixer));
            dat.push_back(color.r);
            dat.push_back(color.g);
            dat.push_back(color.b);
        }
        glBufferData(GL_ARRAY_BUFFER, dat.size() * sizeof(float), dat.data(), GL_DYNAMIC_DRAW);

        glBindVertexArray(VAO);

        // glEnableVertexAttribArray(2);
        // glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)0);
        // glVertexAttribDivisor(2, 1);

        // glEnableVertexAttribArray(3);
        // glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(float)*6, (void*)(sizeof(float)*3));
        // glVertexAttribDivisor(3, 1);

        glDrawElementsInstanced(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (void*) 0, dat.size()/6);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);

        lastTime = newTime;
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}