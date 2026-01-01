#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <algorithm>
#include <string>

// ---------------- CONFIG ----------------
constexpr int N = 256;
constexpr int STEPS = 12000;
constexpr float dt = 1.0f;

// ----------------------------------------

inline int idx(int x, int y) {
    x = (x + N) % N;
    y = (y + N) % N;
    return x + y * N;
}

// Minimal PNG writer (lodepng, public domain)

#include "lodepng.h"

float laplacian(const std::vector<float>& z, int x, int y) {
    return
        -4.0f * z[idx(x,y)]
        + z[idx(x+1,y)]
        + z[idx(x-1,y)]
        + z[idx(x,y+1)]
        + z[idx(x,y-1)];
}

int main() {
    std::mt19937 rng(std::random_device{}());

    std::uniform_real_distribution<float> distDu(0.12f, 0.20f);
    std::uniform_real_distribution<float> distDv(0.04f, 0.10f);
    std::uniform_real_distribution<float> distF (0.015f, 0.040f);
    std::uniform_real_distribution<float> distK (0.045f, 0.070f);

    // Random parameters
    float Du = distDu(rng);
    float Dv = distDv(rng);
    float F  = distF(rng);
    float k  = distK(rng);

    std::cout << "Parameters:\n";
    std::cout << "Du = " << Du << "\n";
    std::cout << "Dv = " << Dv << "\n";
    std::cout << "F  = " << F  << "\n";
    std::cout << "k  = " << k  << "\n";

    // Fields
    std::vector<float> u(N*N, 1.0f);
    std::vector<float> v(N*N, 0.0f);
    std::vector<float> u2(N*N), v2(N*N);

    // Initial perturbation
    for (int y = N/2 - 10; y < N/2 + 10; ++y)
        for (int x = N/2 - 10; x < N/2 + 10; ++x) {
            u[idx(x,y)] = 0.5f;
            v[idx(x,y)] = 0.25f;
        }

    // Simulation loop
    for (int step = 0; step < STEPS; ++step) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                int i = idx(x,y);
                float uvv = u[i] * v[i] * v[i];

                u2[i] = u[i] + dt * (
                    Du * laplacian(u,x,y)
                    - uvv + F * (1 - u[i])
                );

                v2[i] = v[i] + dt * (
                    Dv * laplacian(v,x,y)
                    + uvv - (F + k) * v[i]
                );
            }
        }
        u.swap(u2);
        v.swap(v2);
    }

    // Convert to image (visualize v)
    std::vector<unsigned char> image(N * N * 4);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            float c = std::clamp(v[idx(x,y)], 0.0f, 1.0f);
            unsigned char col = static_cast<unsigned char>(255 * c);

            int p = 4 * (x + y * N);
            image[p+0] = col;
            image[p+1] = col;
            image[p+2] = col;
            image[p+3] = 255;
        }
    }

    std::string filename =
        "turing_Du_" + std::to_string(Du) +
        "_Dv_" + std::to_string(Dv) +
        "_F_"  + std::to_string(F)  +
        "_k_"  + std::to_string(k)  + ".png";

    unsigned error = lodepng::encode(filename, image, N, N);
    if (error)
        std::cerr << "PNG error: " << error << "\n";
    else
        std::cout << "Saved: " << filename << "\n";

    return 0;
}

