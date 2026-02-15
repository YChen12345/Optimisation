#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// Simple support class for a 2D vector
//class vec2D {
//public:
//    float x, y;
//
//    // Default constructor initializes both components to 0
//    vec2D() { x = y = 0.f; };
//
//    // Constructor initializes components with given values
//    vec2D(float _x, float _y) : x(_x), y(_y) {}
//
//    // Constructor initializes components from a vec4
//    vec2D(vec4 v) {
//        x = v[0];
//        y = v[1];
//    }
//
//    // Display the vector components
//    void display() { std::cout << x << '\t' << y << std::endl; }
//
//    // Overloaded subtraction operator for vector subtraction
//    vec2D operator- (vec2D& v) {
//        vec2D q;
//        q.x = x - v.x;
//        q.y = y - v.y;
//        return q;
//    }
//};
class vec2D {
public:
    float x, y;

    vec2D() : x(0), y(0) {}
    vec2D(float _x, float _y) : x(_x), y(_y) {}
    explicit vec2D(const vec4& v) : x(v[0]), y(v[1]) {}

    vec2D operator-(const vec2D& v) const {
        return vec2D(x - v.x, y - v.y);
    }

    float cross(const vec2D& v) const {
        return x * v.y - y * v.x;
    }
};

//Class representing a triangle for rendering purposes
//class triangle {
//    Vertex v[3];       // Vertices of the triangle
//    float area;        // Area of the triangle
//    colour col[3];     // Colors for each vertex of the triangle
//
//public:
//    // Constructor initializes the triangle with three vertices
//    // Input Variables:
//    // - v1, v2, v3: Vertices defining the triangle
//    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
//        v[0] = v1;
//        v[1] = v2;
//        v[2] = v3;
//
//        // Calculate the 2D area of the triangle
//        vec2D e1 = vec2D(v[1].p - v[0].p);
//        vec2D e2 = vec2D(v[2].p - v[0].p);
//        area = std::fabs(e1.x * e2.y - e1.y * e2.x);
//    }
//
//    // Helper function to compute the cross product for barycentric coordinates
//    // Input Variables:
//    // - v1, v2: Edges defining the vector
//    // - p: Point for which coordinates are being calculated
//    float getC(vec2D v1, vec2D v2, vec2D p) {
//        vec2D e = v2 - v1;
//        vec2D q = p - v1;
//        return q.y * e.x - q.x * e.y;
//    }
//
//    // Compute barycentric coordinates for a given point
//    // Input Variables:
//    // - p: Point to check within the triangle
//    // Output Variables:
//    // - alpha, beta, gamma: Barycentric coordinates of the point
//    // Returns true if the point is inside the triangle, false otherwise
//    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
//        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
//        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
//        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;
//
//        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
//        return true;
//    }
//
//    // Template function to interpolate values using barycentric coordinates
//    // Input Variables:
//    // - alpha, beta, gamma: Barycentric coordinates
//    // - a1, a2, a3: Values to interpolate
//    // Returns the interpolated value
//    template <typename T>
//    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
//        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
//    }
//
//    // Draw the triangle on the canvas
//    // Input Variables:
//    // - renderer: Renderer object for drawing
//    // - L: Light object for shading calculations
//    // - ka, kd: Ambient and diffuse lighting coefficients
//    void draw(Renderer& renderer, Light& L, float ka, float kd) {
//        vec2D minV, maxV;
//
//        // Get the screen-space bounds of the triangle
//        getBoundsWindow(renderer.canvas, minV, maxV);
//
//        // Skip very small triangles
//        if (area < 1.f) return;
//
//        // Iterate over the bounding box and check each pixel
//        for (int y = (int)(minV.y); y < (int)ceil(maxV.y); y++) {
//            for (int x = (int)(minV.x); x < (int)ceil(maxV.x); x++) {
//                float alpha, beta, gamma;
//
//                // Check if the pixel lies inside the triangle
//                if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
//                    // Interpolate color, depth, and normals
//                    colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
//                    c.clampColour();
//                    float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
//                    vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
//                    normal.normalise();
//
//                    // Perform Z-buffer test and apply shading
//                    if (renderer.zbuffer(x, y) > depth && depth > 0.001f) {
//                        // typical shader begin
//                        L.omega_i.normalise();
//                        float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
//                        colour a = (c * kd) * (L.L * dot) + (L.ambient * ka); // using kd instead of ka for ambient
//                        // typical shader end
//                        unsigned char r, g, b;
//                        a.toRGB(r, g, b);
//                        renderer.canvas.draw(x, y, r, g, b);
//                        renderer.zbuffer(x, y) = depth;
//                    }
//                }
//            }
//        }
//    }
//
//    // Compute the 2D bounds of the triangle
//    // Output Variables:
//    // - minV, maxV: Minimum and maximum bounds in 2D space
//    void getBounds(vec2D& minV, vec2D& maxV) {
//        minV = vec2D(v[0].p);
//        maxV = vec2D(v[0].p);
//        for (unsigned int i = 1; i < 3; i++) {
//            minV.x = std::min(minV.x, v[i].p[0]);
//            minV.y = std::min(minV.y, v[i].p[1]);
//            maxV.x = std::max(maxV.x, v[i].p[0]);
//            maxV.y = std::max(maxV.y, v[i].p[1]);
//        }
//    }
//
//    // Compute the 2D bounds of the triangle, clipped to the canvas
//    // Input Variables:
//    // - canvas: Reference to the rendering canvas
//    // Output Variables:
//    // - minV, maxV: Clipped minimum and maximum bounds
//    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
//        getBounds(minV, maxV);
//        minV.x = std::max(minV.x, static_cast<float>(0));
//        minV.y = std::max(minV.y, static_cast<float>(0));
//        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
//        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
//    }
//
//    // Debugging utility to display the triangle bounds on the canvas
//    // Input Variables:
//    // - canvas: Reference to the rendering canvas
//    void drawBounds(GamesEngineeringBase::Window& canvas) {
//        vec2D minV, maxV;
//        getBounds(minV, maxV);
//
//        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
//            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
//                canvas.draw(x, y, 255, 0, 0);
//            }
//        }
//    }
//
//    // Debugging utility to display the coordinates of the triangle vertices
//    void display() {
//        for (unsigned int i = 0; i < 3; i++) {
//            v[i].p.display();
//        }
//        std::cout << std::endl;
//    }
//};

class triangle {
    Vertex v[3];

    int minX = 0, minY = 0, maxX = -1, maxY = -1;

    float A0 = 0, B0 = 0, C0 = 0;
    float A1 = 0, B1 = 0, C1 = 0;
    float A2 = 0, B2 = 0, C2 = 0;

    float invArea = 0.0f;

    static inline void edgeCoeffs(const Vertex& a, const Vertex& b, float& A, float& B, float& C) noexcept {
        A = (a.p[1] - b.p[1]);
        B = (b.p[0] - a.p[0]);
        C = (a.p[0] * b.p[1] - b.p[0] * a.p[1]);
    }

    static inline float evalEdge(float A, float B, float C, float x, float y) noexcept {
        return A * x + B * y + C;
    }

    static inline float area2(const Vertex& a, const Vertex& b, const Vertex& c) noexcept {
        return (b.p[0] - a.p[0]) * (c.p[1] - a.p[1]) - (b.p[1] - a.p[1]) * (c.p[0] - a.p[0]);
    }

    inline void buildBounds(const GamesEngineeringBase::Window& canvas) noexcept {
        float fx0 = v[0].p[0], fy0 = v[0].p[1];
        float fx1 = v[1].p[0], fy1 = v[1].p[1];
        float fx2 = v[2].p[0], fy2 = v[2].p[1];

        float fminX = std::min({ fx0, fx1, fx2 });
        float fminY = std::min({ fy0, fy1, fy2 });
        float fmaxX = std::max({ fx0, fx1, fx2 });
        float fmaxY = std::max({ fy0, fy1, fy2 });

        int W = canvas.getWidth();
        int H = canvas.getHeight();

        minX = std::max(0, (int)std::floor(fminX));
        minY = std::max(0, (int)std::floor(fminY));
        maxX = std::min(W - 1, (int)std::ceil(fmaxX));
        maxY = std::min(H - 1, (int)std::ceil(fmaxY));
    }

public:
    triangle(const Vertex& v0, const Vertex& v1, const Vertex& v2) {
        v[0] = v0; v[1] = v1; v[2] = v2;
    }

    void draw(Renderer& renderer, const Light& Lin, float ka, float kd) {
        drawBand(renderer, Lin, ka, kd, 0, renderer.canvas.getHeight());
    }

    void drawBand(Renderer& renderer, const Light& Lin, float ka, float kd, int y0, int y1) {
        buildBounds(renderer.canvas);
        if (maxX < minX || maxY < minY) return;

        float a2 = area2(v[0], v[1], v[2]);
        if (a2 <= 1e-6f) return;
        invArea = 1.0f / a2;

        edgeCoeffs(v[1], v[2], A0, B0, C0);
        edgeCoeffs(v[2], v[0], A1, B1, C1);
        edgeCoeffs(v[0], v[1], A2, B2, C2);

        int startY = std::max(minY, y0);
        int endY = std::min(maxY + 1, y1);

        if (startY >= endY) return;
        vec4 Ldir = Lin.omega_i;
        Ldir.normalise();

        for (int y = startY; y < endY; ++y) {
            float py = (float)y + 0.5f;
            float px0 = (float)minX + 0.5f;

            float w0 = evalEdge(A0, B0, C0, px0, py);
            float w1 = evalEdge(A1, B1, C1, px0, py);
            float w2 = evalEdge(A2, B2, C2, px0, py);

            for (int x = minX; x <= maxX; ++x) {
                if (w0 >= 0.0f && w1 >= 0.0f && w2 >= 0.0f) {
                    float b0 = w0 * invArea;
                    float b1 = w1 * invArea;
                    float b2 = w2 * invArea;

                    float depth = v[0].p[2] * b0 + v[1].p[2] * b1 + v[2].p[2] * b2;
                    if (depth > 0.001f && renderer.zbuffer(x, y) > depth) {
                        vec4 n = v[0].normal * b0 + v[1].normal * b1 + v[2].normal * b2;
                        n.normalise();

                        float ndotl = std::max(vec4::dot(Ldir, n), 0.0f);

                        colour c = v[0].rgb * b0 + v[1].rgb * b1 + v[2].rgb * b2;

                        float diffuse = kd * ndotl;

                        colour out = (c * diffuse) + (Lin.ambient*ka);

                        unsigned char r, g, b;
                        out.toRGB(r, g, b);

                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                    }
                }

                w0 += A0;
                w1 += A1;
                w2 += A2;
            }
        }
    }
};