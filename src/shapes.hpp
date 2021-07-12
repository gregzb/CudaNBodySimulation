#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>

template <class T>
void generateSphereMesh(std::vector<T> &vertices, std::vector<unsigned int> &indices, int thetaSteps, int phiSteps)
{
    // radius = 1, origin = (0, 0, 0)

    double thetaStepSize = 2 * M_PI / thetaSteps;
    double phiStepSize = M_PI / phiSteps;

    std::vector<glm::vec3> points((phiSteps-1) * thetaSteps + 2, {0, 0, 0});

    auto idx_of = [phiSteps, thetaSteps](int phiStep, int thetaStep) {
        thetaStep %= thetaSteps;
        if (phiStep == 0) return 0;
        // if (phiStep == phiSteps) return (phiSteps-1-1) * thetaSteps + (thetaSteps-1) + 1 + 1;
        if (phiStep == phiSteps) return (phiSteps-1) * thetaSteps + 1;
        return (phiStep-1) * thetaSteps + thetaStep + 1;
    };

    // std::vector<std::vector<glm::vec3>> dp(phiSteps+1, std::vector<glm::vec3>(thetaSteps, {0, 0, 0}));
    for (int phiStep = 0; phiStep <= phiSteps; phiStep++)
    {
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++)
        {
            double phi = phiStep * phiStepSize;
            double theta = thetaStep * thetaStepSize;
            // dp[phiStep][thetaStep] = {std::cos(phi),
            //                           std::sin(phi) * std::cos(theta),
            //                           std::sin(phi) * std::sin(theta)};
            int idx = idx_of(phiStep, thetaStep);
            points[idx] = {std::cos(phi),
                           std::sin(phi) * std::cos(theta),
                           std::sin(phi) * std::sin(theta)};
        }
    }

    vertices.clear();
    for (unsigned i = 0; i < points.size(); i++) {
        vertices.emplace_back(points[i], points[i]);
    }

    //merge all thetasteps into one at phistep = 0 and phistep = phiSteps

    indices.clear();

    const std::vector<int> offsets{0, 0, 1, 0, 1, 1, 0, 1};

    for (int phiStep = 0; phiStep < phiSteps; phiStep++)
    {
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++)
        {
            if (phiStep != 0) {
                indices.push_back(idx_of(phiStep, thetaStep));
                indices.push_back(idx_of(phiStep+1, thetaStep+1));
                indices.push_back(idx_of(phiStep, thetaStep+1));
            }
            if (phiStep != phiSteps-1) {
                indices.push_back(idx_of(phiStep, thetaStep));
                indices.push_back(idx_of(phiStep+1, thetaStep));
                indices.push_back(idx_of(phiStep+1, thetaStep+1));
            }
        }
    }
}