#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>
#include <thread>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"
#include <vector>

struct FinalTri {
    triangle tri;
    float ka;
    float kd;
    int minY;
    int maxY; // [minY, maxY] triangle Y coverage in screen space

    FinalTri(const triangle& t, float _ka, float _kd, int _minY, int _maxY)
        : tri(t), ka(_ka), kd(_kd), minY(_minY), maxY(_maxY) {
    }
};

void buildTris(Renderer& renderer, Mesh* mesh, matrix& camera, std::vector<FinalTri>& out) {
    matrix p = renderer.perspective * camera * mesh->world;

    const int H = renderer.canvas.getHeight();

    for (triIndices& ind : mesh->triangles) {
        Vertex t[3];

        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p;
            t[i].p.divideW();

            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
            t[i].normal.normalise();

            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1];

            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        float fy0 = std::min({ t[0].p[1], t[1].p[1], t[2].p[1] });
        float fy1 = std::max({ t[0].p[1], t[1].p[1], t[2].p[1] });
        int minY = (int)std::floor(fy0);
        int maxY = (int)std::ceil(fy1);

        if (maxY < 0 || minY >= H) continue;
        if (minY < 0) minY = 0;
        if (maxY >= H) maxY = H - 1;

        out.emplace_back(triangle(t[0], t[1], t[2]), mesh->ka, mesh->kd, minY, maxY);
    }
}

int chooseThreadCount(int maxThreads, int H, size_t jobCount) 
{
    int hw = (int)std::thread::hardware_concurrency();
    if (hw <= 0) hw = maxThreads;

    int n = std::min(maxThreads, std::min(hw, H));

    if (jobCount < 1500) n = std::min(n, 2);
    else if (jobCount < 4000) n = std::min(n, 4);

    return std::max(1, n);
}

static void renderMT(Renderer& renderer,const std::vector<Mesh*>& meshes,matrix& camera,Light& L,int maxThreads = 11) {
    size_t totalTris = 0;
    for (auto* m : meshes) totalTris += m->triangles.size();

    std::vector<FinalTri> tris;
    tris.reserve(totalTris);

    for (auto* m : meshes) buildTris(renderer, m, camera, tris);

    const int H = renderer.canvas.getHeight();
    const int nThreads = chooseThreadCount(maxThreads, H, tris.size());
    const int band = (H + nThreads - 1) / nThreads;

    std::vector<std::thread> workers;
    workers.reserve(nThreads);

    for (int t = 0; t < nThreads; ++t) {
        const int y0 = t * band;
        const int y1 = std::min(H, y0 + band);
        if (y0 >= y1) break;

        workers.emplace_back([&renderer, &L, &tris, y0, y1]() {
            for (auto& j : tris) {
                if (j.maxY < y0 || j.minY >= y1) continue;
                j.tri.drawBand(renderer, L, j.ka, j.kd, y0, y1);
            }
            });
    }

    for (auto& th : workers) th.join();
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------
static int scene_state = 0;
// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    int roundCount = 0;
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    bool running = true;

    std::vector<Mesh*> scene;

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    float timer = 0;
    int frameCount = 0;

    // Main rendering loop
    while (running) {
        auto frame_start = std::chrono::high_resolution_clock::now();
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << "scene1:" << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                std::cout << "scene1_deltaTime_average" << ":" <<timer/frameCount<< "ms\n";
                std::cout << "scene1_FPS_average" << ":" << 1000*frameCount/timer << "fps\n";
                start = std::chrono::high_resolution_clock::now();
                roundCount++;
                if (roundCount > 1) {
                    std::cout << "\n";
                    scene_state++;
                    break;
                }
            }
        }

        /* for (auto& m : scene) {
            render(renderer, m, camera, L);
        } */     
        renderMT(renderer, scene, camera, L, 2);
        
        renderer.present();
        auto frame_end = std::chrono::high_resolution_clock::now();
        timer += std::chrono::duration<double, std::milli>(frame_end - frame_start).count();
        frameCount++;
        //std::cout << "scene1_deltaTime" << ":" << std::chrono::duration<double, std::milli>(frame_end - frame_start).count() << "ms\n";
        //std::cout << "scene1_FPS" << ":" << static_cast<int>(1000 / std::chrono::duration<double, std::milli>(frame_end - frame_start).count()) << "fps\n";
    }

    for (auto& m : scene)
        delete m;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------
// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    int roundCount = 0;
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.push_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    float timer = 0;
    int frameCount = 0;

    bool running = true;
    while (running) {
        auto frame_start = std::chrono::high_resolution_clock::now();
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << "scene2:" << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                std::cout << "scene2_deltaTime_average" << ":" << timer / frameCount << "ms\n";
                std::cout << "scene2_FPS_average" << ":" << 1000 * frameCount / timer << "fps\n";
                start = std::chrono::high_resolution_clock::now();
                roundCount++;
                if (roundCount > 1) {
                    std::cout <<"\n";
                    scene_state++;
                    break;
                }
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        /*for (auto& m : scene) {
            render(renderer, m, camera, L);
        }*/
        renderMT(renderer, scene, camera, L, 6);

        renderer.present();
        auto frame_end = std::chrono::high_resolution_clock::now();
        timer += std::chrono::duration<double, std::milli>(frame_end - frame_start).count();
        frameCount++;
        //std::cout << "scene2_deltaTime" << ":" << std::chrono::duration<double, std::milli>(frame_end - frame_start).count() << "ms\n";
        //std::cout << "scene2_FPS" << ":" << static_cast<int>(1000 / std::chrono::duration<double, std::milli>(frame_end - frame_start).count()) << "fps\n";
    }

    for (auto& m : scene) {
        delete m;
    }          
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    int roundCount = 0;
    Renderer renderer;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };
    matrix camera = matrix::makeIdentity();
    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; };
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    for (unsigned int z = 0; z < 10; z++) {
        for (unsigned int y = 0; y < 6; y++) {
            for (unsigned int x = 0; x < 8; x++) {
                Mesh* m = new Mesh();
                *m = Mesh::makeCube(1.f);
                scene.push_back(m);
                m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -6.f*z);
                rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
                rotations.push_back(r);
            }
        }
    }    
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        //--------------------------------------------------------------------------------------------------------------
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        //--------------------------------------------------------------------------------------------------------------
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(0.0f, 2.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        //---------------------------------------------------------------------------------------------------------------
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(0.0f, -2.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    bool running = true;

    //------------------------------------------------------------
    Mesh* mesh = new Mesh();
    *mesh = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(mesh);
    //------------------------------------------------------------
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(0.f, sphereOffset, -30.f);
    //------------------------------------------------------------------------------
    float x = 0.0f, y = 0.0f, z = -4.0f;
    mesh->world = matrix::makeTranslation(x, y, z);
    //------------------------------------------------------------------------------

    float zoffset = 8.0f;
    float step = -0.1f;
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    float timer = 0;
    int frameCount = 0;

    while (running) {
        auto frame_start = std::chrono::high_resolution_clock::now();
        renderer.canvas.checkInput();
        renderer.clear();
        zoffset += step;
        sphereOffset += sphereStep;
        //-----------------------------------------------------
        camera = matrix::makeTranslation(0, 0, -zoffset);
        //-----------------------------------------------------
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
        mesh->world = matrix::makeTranslation(x, y, zoffset+z);
        sphere->world = matrix::makeTranslation(0.f, sphereOffset, -30.f);
        //----------------------------------------------------------------
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << "sceneTest:" << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                std::cout << "sceneTest_deltaTime_average" << ":" << timer / frameCount << "ms\n";
                std::cout << "sceneTest_FPS_average" << ":" << 1000 * frameCount / timer << "fps\n";
                start = std::chrono::high_resolution_clock::now();
                roundCount++;
                if (roundCount > 1) {
                    std::cout << "\n";
                    scene_state++;
                    break;
                }
            }
        }
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
        }
        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;
        if (x > 4) x = 4;
        if (x < -4)  x = -4;
        if (y > 4) y = 4;
        if (y < -4)  y = -4;
        if (z > -4) z = -4;
        if (z < -10) z = -10;

        // Render each object in the scene
        /*for (auto& m : scene) {
            render(renderer, m, camera, L);
        }*/
        renderMT(renderer, scene, camera, L, 11);

        renderer.present(); // Display the rendered frame
        auto frame_end = std::chrono::high_resolution_clock::now();
        timer += std::chrono::duration<double, std::milli>(frame_end - frame_start).count();
        frameCount++;
        //std::cout << "sceneTest_deltaTime" << ":" << std::chrono::duration<double, std::milli>(frame_end - frame_start).count() << "ms\n";
        //std::cout << "sceneTest_FPS" << ":" << static_cast<int>(1000 / std::chrono::duration<double, std::milli>(frame_end - frame_start).count()) << "fps\n";
    }
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------
// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    if (scene_state == 0) {
        scene1();
    }
    if (scene_state == 1) {
        scene2();
    }
    if (scene_state == 2) {
        sceneTest();
    }
    return 0;
}