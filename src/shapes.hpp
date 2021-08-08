#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>

// template <class T>
void generateSphereMesh(std::vector<glm::vec3> &vertices, std::vector<unsigned int> &indices, int thetaSteps, int phiSteps)
{
    // radius = 1, origin = (0, 0, 0)

    double thetaStepSize = 2 * M_PI / thetaSteps;
    double phiStepSize = M_PI / phiSteps;

    std::vector<glm::vec3> points((phiSteps-1) * thetaSteps + 2, {0, 0, 0});

    auto idx_of = [phiSteps, thetaSteps](int phiStep, int thetaStep) {
        thetaStep %= thetaSteps;
        if (phiStep == 0) return 0;
        if (phiStep == phiSteps) return (phiSteps-1) * thetaSteps + 1;
        return (phiStep-1) * thetaSteps + thetaStep + 1;
    };

    for (int phiStep = 0; phiStep <= phiSteps; phiStep++)
    {
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++)
        {
            double phi = phiStep * phiStepSize;
            double theta = thetaStep * thetaStepSize;
            int idx = idx_of(phiStep, thetaStep);
            points[idx] = {std::cos(phi),
                           std::sin(phi) * std::cos(theta),
                           std::sin(phi) * std::sin(theta)};
        }
    }

    vertices.clear();
    for (unsigned i = 0; i < points.size(); i++) {
        vertices.emplace_back(points[i]);
    }

    //merge all thetasteps into one at phistep = 0 and phistep = phiSteps

    indices.clear();
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

void generateTorusMesh(std::vector<glm::vec3> &vertices, std::vector<unsigned int> &indices, double r1, double r2, int thetaSteps, int phiSteps)
{
    double thetaStepSize = 2 * M_PI / thetaSteps;
    double phiStepSize = 2 * M_PI / phiSteps;

    std::vector<glm::vec3> points(thetaSteps*phiSteps);

    auto idx_of = [phiSteps, thetaSteps](int phiStep, int thetaStep) {
        return phiStep * thetaSteps + thetaStep;
    };

    for (int phiStep = 0; phiStep < phiSteps; phiStep++)
    {
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++)
        {
            double phi = phiStep * phiStepSize;
            double theta = thetaStep * thetaStepSize;
            points[idx_of(phiStep, thetaStep)] = {(r2 + r1 * std::cos(phi)) * std::cos(theta),
                            r1 * std::sin(phi),
                            (r2 + r1 * std::cos(phi)) * std::sin(theta)};
        }
    }

    vertices.clear();
    for (unsigned i = 0; i < points.size(); i++) {
        vertices.emplace_back(points[i]);
    }

    indices.clear();

    for (int phiStep = 0; phiStep < phiSteps; phiStep++)
    {
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++)
        {
                indices.push_back(idx_of(phiStep, thetaStep));
                indices.push_back(idx_of((phiStep+1)%phiSteps, (thetaStep+1)%thetaSteps));
                indices.push_back(idx_of(phiStep, (thetaStep+1)%thetaSteps));
                
                indices.push_back(idx_of(phiStep, thetaStep));
                indices.push_back(idx_of(phiStep, (thetaStep+1)%thetaSteps));
                indices.push_back(idx_of((phiStep+1)%phiSteps, thetaStep));
        }
    }
}