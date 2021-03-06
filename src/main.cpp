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

#include <filesystem>
#include <fenv.h>

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void process_input(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
    screen_width = width;
    screen_height = height;
}

GLFWwindow *initialize_glfw() {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_FALSE); // single buffering, inherently disables vsync?

    // glfw window creation
    // --------------------
    GLFWwindow *window = glfwCreateWindow(screen_width, screen_height, "N-Body Simulation", NULL, NULL);
    if (window == NULL) {
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    return window;
}

void init_sphere_buffers(unsigned int &VBO, unsigned int &VAO, unsigned int &EBO, const std::vector<vertex> &vertices,
                         const std::vector<unsigned int> &indices) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void *) 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void *) (sizeof(float) * 3));
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
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
}

// void init_bodies(std::vector<body> &bodies) {
//     // std::random_device rd;
//     // std::mt19937 gen(rd());
//     std::mt19937 gen(0); // keep it reproducible for now
//     std::uniform_real_distribution<float> radius_distrib(1, 2);
//     std::uniform_real_distribution<float> phi_distrib(0, (float)M_PI);
//     std::uniform_real_distribution<float> theta_distrib(0, (float)M_PI*2);

//     glm::vec3 center(0, 0, 0);
//     bodies.push_back({center, {0, 0, 0}, 10000000});

//     float initial_velocity_mag = 0.02;

//     for (int i = 0; i < 100; i++) {
//         float r = radius_distrib(gen);
//         float phi = phi_distrib(gen);
//         float theta = theta_distrib(gen);
//         glm::vec3 pos(r*std::cos(theta)*std::sin(phi), r*std::sin(theta)*std::sin(phi), r*std::cos(phi));
//         glm::vec3 normal(glm::normalize(pos-center));

//         float epsilon = 0.01;

//         glm::vec3 another;
//         while(glm::length2(cross(another=glm::sphericalRand(1.0f), normal)) < epsilon);

//         glm::vec3 vector_on_plane = glm::normalize(cross(normal, another)) * initial_velocity_mag;

//         bodies.push_back({pos, vector_on_plane, 0.01});
//     }
// }

//void init_bodies(std::vector<body> &bodies) {
//    std::mt19937 gen(0); // keep it reproducible for now
//    std::uniform_real_distribution<float> radius_distrib(1, 2);
//    std::uniform_real_distribution<float> phi_distrib(0, (float)M_PI);
//    std::uniform_real_distribution<float> theta_distrib(0, (float)M_PI*2);
//
//    std::uniform_real_distribution<float> square_distrib(-0.4, 0.4);
//
//    glm::vec3 center(0, 0, 0);
//    bodies.push_back({center, {0, 0, 0}, 10000000});
//
//    float initial_velocity_mag = 0.02;
//
//    for (int i = 0; i < 200000; i++) {
//        bodies.push_back({{square_distrib(gen)+1.4, square_distrib(gen), square_distrib(gen)}, {0, initial_velocity_mag, 0}, 0.01f});
//    }
//}

void add_circular_body_system(std::vector<body> &bodies, const glm::vec3 &center, const glm::vec3 &center_velocity,
                              int num_orbiting) {
    std::mt19937 gen(0); // keep it reproducible for now
    std::uniform_real_distribution<float> radius_distrib(0.4, 1.8);
    std::uniform_real_distribution<float> theta_distrib(0, (float) M_PI * 2);

    std::uniform_real_distribution<float> square_distrib(-0.1, 0.1);

    bodies.push_back({center, center_velocity, 10000000});

    float initial_velocity_mag_multiplier = 0.02;

    for (int i = 0; i < num_orbiting; i++) {
        float r = radius_distrib(gen);
        float theta = theta_distrib(gen);

        float x_pos = r * std::cos(theta);
        float y_pos = r * std::sin(theta);
        float z_pos = square_distrib(gen);

        glm::vec3 pos{x_pos, y_pos, z_pos};
        glm::vec3 velocity{-y_pos, x_pos, 0};
        velocity = glm::normalize(velocity);
        velocity *= initial_velocity_mag_multiplier * std::pow(2.6 - r, 0.6f);

        bodies.push_back({center + pos, velocity, 0.01f});
    }
}

void init_bodies(std::vector<body> &bodies) {
    add_circular_body_system(bodies, {0, -4, -0.1}, {0.006, 0, 0}, 75000);
    add_circular_body_system(bodies, {0, 4, 0.1}, {-0.006, 0, 0}, 75000);
}

void imgui_setup(GLFWwindow *window) {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void) io;

    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");
}

shader_program make_shader_program(const std::string &vertex_shader_file, const std::string &fragment_shader_file) {
    shader<GL_VERTEX_SHADER> vertex_shader(vertex_shader_file);
    shader<GL_FRAGMENT_SHADER> fragment_shader(fragment_shader_file);
    shader_program shader_prog({vertex_shader.get_id(), fragment_shader.get_id()});
    vertex_shader.release();
    fragment_shader.release();

    return shader_prog;
}

int main() {
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    std::cout << "Starting application with cwd: " << std::filesystem::current_path() << std::endl;

    glfwInit();

    GLFWwindow *window;
    if (!(window = initialize_glfw())) {
        std::cout << "Failed to create GLFW window" << std::endl;
        return -1;
    }

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    imgui_setup(window);

    shader_program body_shader = make_shader_program("shaders/body_vertex.glsl", "shaders/body_fragment.glsl");
    body_shader.add_uniform("model");
    body_shader.add_uniform("view");
    body_shader.add_uniform("projection");
    body_shader.add_uniform("normal_model_view");
    body_shader.add_uniform("light_pos");

    shader_program highlighter_shader = make_shader_program("shaders/highlighter_vertex.glsl",
                                                            "shaders/highlighter_fragment.glsl");
    highlighter_shader.add_uniform("screen_size");
    highlighter_shader.add_uniform("highlight_center");
    highlighter_shader.add_uniform("radii");
    highlighter_shader.add_uniform("highlight_color");

    std::vector<vertex> vertices;
    std::vector<unsigned int> indices;

    std::vector<glm::vec3> vertex_locations;
    generateSphereMesh(vertex_locations, indices, 8, 8);

    for (unsigned i = 0; i < vertex_locations.size(); i++) {
        // position, normal (they're the same for a sphere)
        vertices.emplace_back(vertex_locations[i], vertex_locations[i]);
    }

    unsigned int VBO, VAO, EBO;
    init_sphere_buffers(VBO, VAO, EBO, vertices, indices);

    std::vector<body> bodies;
    init_bodies(bodies);

    // std::vector<glm::vec2> screen_vertices{{-1, -1}, {-1, 1}, {1, 1}, {1, -1}};
    std::vector<float> screen_vertices{-1, -1, 1, -1, 1, 1, -1, 1};
    std::vector<unsigned int> screen_indices{0, 1, 2, 0, 2, 3};

    unsigned int screen_VBO, screen_VAO, screen_EBO;
    glGenVertexArrays(1, &screen_VAO);
    glGenBuffers(1, &screen_VBO);
    glGenBuffers(1, &screen_EBO);

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(screen_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, screen_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * screen_vertices.size(), screen_vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screen_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * screen_indices.size(), screen_indices.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, (void *) 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);


    unsigned int instance_buffer;
    glGenBuffers(1, &instance_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, instance_buffer);
    std::vector<float> dat;
    for (const auto &body : bodies) {
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
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void *) 0);
    glVertexAttribDivisor(2, 1);

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void *) (sizeof(float) * 3));
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(0);

    enable_gl_settings();

    float time_scale = 0.1f;
    nbody_simulation simulation(bodies, time_scale);

    camera cam(glm::radians(60.0f), {0, 0, 0}, {0, 0, 0});

    double lastTime = glfwGetTime();

    double frame_time = 0;

    // std::random_device rd; // obtain a random number from hardware
    // std::mt19937 gen(rd()); // seed the generator
    std::mt19937 gen(0);
    std::uniform_int_distribution<> distr(0, bodies.size() - 1); // define the range

    int selected_body_idx = -1;

    ImVec4 highlight_color(0.2f, 0.9f, 0.2f, 0.5f);
    glm::vec2 highlight_radii(0.01f, 0.02f);

    float body_size = 1.0f;

    int backend_int = static_cast<int>(nbody_simulation::CalculationBackend::NAIVE_CPU);
    float barnes_hut_factor = 1.0f;

    bool sim_running = false;
    bool prev_sim_running = sim_running;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window)) {
        double newTime = glfwGetTime();
        double dt = newTime - lastTime;
        while (dt < (1 / 90.0)) {
            newTime = glfwGetTime();
            dt = newTime - lastTime;
        }

        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            ImGui::Begin("Simulation Control Panel");

            if (ImGui::CollapsingHeader("Simulation Settings")) {
                ImGui::SliderFloat("Time Step", &time_scale, 0.001f, 100.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::SliderFloat("Body Size", &body_size, 0.1f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);

                ImGui::Dummy(ImVec2(0.0f, 5.0f));

                ImGui::RadioButton("CPU - Naive", &backend_int,
                                   static_cast<int>(nbody_simulation::CalculationBackend::NAIVE_CPU));
                ImGui::SameLine();
                ImGui::RadioButton("CPU - Barnes Hut", &backend_int,
                                   static_cast<int>(nbody_simulation::CalculationBackend::BARNES_HUT_CPU));
                ImGui::RadioButton("GPU - Naive", &backend_int,
                                   static_cast<int>(nbody_simulation::CalculationBackend::NAIVE_GPU));
                ImGui::SameLine();
                ImGui::RadioButton("GPU - Barnes Hut", &backend_int,
                                   static_cast<int>(nbody_simulation::CalculationBackend::BARNES_HUT_GPU));
                ImGui::SliderFloat("Barnes Hut Factor", &barnes_hut_factor, 0.1f, 5.0f, "%.2f");
            }

            ImGui::Dummy(ImVec2(0.0f, 5.0f));

            if (ImGui::CollapsingHeader("Highlight Settings")) {
                ImGui::ColorEdit4("Highlight Color", (float *) &highlight_color);
                ImGui::SliderFloat("Inner Radius", &highlight_radii.x, 0.0f, 0.25f, "%.3f",
                                   ImGuiSliderFlags_Logarithmic);
                ImGui::SliderFloat("Outer Radius", &highlight_radii.y, 0.0f, 0.25f, "%.3f",
                                   ImGuiSliderFlags_Logarithmic);

                if (ImGui::Button("Highlight Random Body")) {
                    selected_body_idx = distr(gen);
                }
                ImGui::SameLine();
                if (ImGui::Button("Unselect Body")) {
                    selected_body_idx = -1;
                }
            }

            ImGui::Dummy(ImVec2(0.0f, 5.0f));

            std::string sim_action = sim_running ? "Stop Simulation" : "Launch Simulation";

            if (ImGui::Button(sim_action.c_str())) {
                sim_running = !sim_running;
            }

            ImGui::Text("Application averaging %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
                        ImGui::GetIO().Framerate);
            ImGui::End();
        }

        process_input(window);

        // render
        // ------
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        body_shader.use();

        cam.set_pos({0, 0, 9});
//        cam.set_pos({glm::sin(newTime) * 5.5f, 0, glm::cos(newTime)*5.5f});
//        cam.set_target({0, 0, 0});

        glm::mat4 model(1.0f);
        float base_size = 0.03f;
        float scaled_size = base_size * body_size;
        model = glm::scale(model, glm::vec3(scaled_size));
        const glm::mat4 &view = cam.get_view_matrix();
        const glm::mat4 &projection = cam.get_projection_matrix();

        body_shader.set_mvp(model, view, projection);
        glm::mat3 normal_model_view = glm::transpose(glm::inverse(glm::mat3(cam.get_view_matrix() * model)));
        glUniformMatrix3fv(body_shader.get_uniform_location("normal_model_view"), 1, GL_FALSE,
                           &normal_model_view[0][0]);

        glm::vec3 light_world_pos{0, 8, 8};
        glm::vec3 light_pos = glm::vec3(view * glm::vec4(light_world_pos, 1.0f));
        glUniform3f(body_shader.get_uniform_location("light_pos"), light_pos.x, light_pos.y, light_pos.z);

        glBindVertexArray(VAO);

        simulation.set_time_scale(time_scale);
        simulation.set_barnes_hut_factor(barnes_hut_factor);
        simulation.set_backend(static_cast<nbody_simulation::CalculationBackend>(backend_int));

        if (sim_running) {
            if (!prev_sim_running) {
                simulation.init();
            }
            simulation.step();
        }

        glm::vec3 red(1, 0.2, 0.2);
        glm::vec3 blue(0.2, 0.4, 1);

        glBindBuffer(GL_ARRAY_BUFFER, instance_buffer);
        dat.clear();
        int num_bodies = simulation.get_bodies().size();
        for (int i = 0; i < num_bodies; i++) {
            const auto &body = simulation.get_bodies()[i];
            dat.push_back(body.pos.x);
            dat.push_back(body.pos.y);
            dat.push_back(body.pos.z);
            float mag = glm::length(body.vel);
            float mixer = glm::clamp(mag * 25, 0.0f, 1.0f);
            glm::vec3 color(glm::mix(blue, red, mixer));
            if (i == selected_body_idx) {
                dat.push_back(highlight_color.x * 1.25f);
                dat.push_back(highlight_color.y * 1.25f);
                dat.push_back(highlight_color.z * 1.25f);
            } else {
                dat.push_back(color.r);
                dat.push_back(color.g);
                dat.push_back(color.b);
            }
        }
        glBufferData(GL_ARRAY_BUFFER, dat.size() * sizeof(float), dat.data(), GL_DYNAMIC_DRAW);

        glDrawElementsInstanced(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (void *) 0, dat.size() / 6);

        glm::vec4 target_location{-2, -2, 0, 1};
        if (selected_body_idx >= 0) {
            target_location = projection * view * glm::vec4(simulation.get_bodies()[selected_body_idx].pos, 1.0f);
            target_location /= target_location.w;
        }

        highlighter_shader.use();
        glUniform2f(highlighter_shader.get_uniform_location("screen_size"), screen_width, screen_height);
        glUniform2f(highlighter_shader.get_uniform_location("highlight_center"), target_location.x, target_location.y);
        glUniform2f(highlighter_shader.get_uniform_location("radii"), highlight_radii.x, highlight_radii.y);
        glUniform4f(highlighter_shader.get_uniform_location("highlight_color"), highlight_color.x, highlight_color.y,
                    highlight_color.z, highlight_color.w);
        glBindVertexArray(screen_VAO);
        glClear(GL_DEPTH_BUFFER_BIT);
        glDrawElements(GL_TRIANGLES, screen_indices.size(), GL_UNSIGNED_INT, (void *) 0);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);

        lastTime = newTime;
        prev_sim_running = sim_running;
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}