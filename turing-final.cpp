#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <ctime>

class TuringPattern {
private:
    int width, height;
    double* u;
    double* v;
    double* u_next;
    double* v_next;
    
    double Du, Dv;
    double feed, kill;
    double dt;
    
    unsigned int numThreads;
    std::string filename;
    
    inline int idx(int y, int x) const {
        return y * width + x;
    }
    
    std::string generateFilename() {
        std::ostringstream oss;
        oss << "turing_"
            << "Du" << std::fixed << std::setprecision(4) << Du << "_"
            << "Dv" << std::setprecision(4) << Dv << "_"
            << "f" << std::setprecision(4) << feed << "_"
            << "k" << std::setprecision(4) << kill << "_"
            << std::time(nullptr);
        return oss.str();
    }
    
public:
    TuringPattern(int w, int h, double du, double dv, double f, double k) 
        : width(w), height(h), Du(du), Dv(dv), feed(f), kill(k), dt(1.0) {
        
        numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
        
        filename = generateFilename();
        
        size_t size = width * height;
        u = new double[size];
        v = new double[size];
        u_next = new double[size];
        v_next = new double[size];
        
        // Initialize arrays
        for (size_t i = 0; i < size; i++) {
            u[i] = 1.0;
            v[i] = 0.0;
        }
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        int cx = width / 2;
        int cy = height / 2;
        int radius = std::min(width, height) / 10;
        
        for (int i = -radius; i <= radius; i++) {
            for (int j = -radius; j <= radius; j++) {
                if (i*i + j*j <= radius*radius) {
                    int y = cy + i;
                    int x = cx + j;
                    if (y >= 0 && y < height && x >= 0 && x < width) {
                        int index = idx(y, x);
                        u[index] = 0.5 + 0.1 * dis(gen);
                        v[index] = 0.25 + 0.1 * dis(gen);
                    }
                }
            }
        }
    }
    
    ~TuringPattern() {
        delete[] u;
        delete[] v;
        delete[] u_next;
        delete[] v_next;
    }
    
    inline double laplacian(const double* grid, int y, int x) const {
        int ym = (y - 1 + height) % height;
        int yp = (y + 1) % height;
        int xm = (x - 1 + width) % width;
        int xp = (x + 1) % width;
        
        return grid[idx(y, xp)] + grid[idx(y, xm)] + 
               grid[idx(yp, x)] + grid[idx(ym, x)] - 
               4.0 * grid[idx(y, x)];
    }
    
    void stepRange(int startY, int endY) {
        const double dt_Du = dt * Du;
        const double dt_Dv = dt * Dv;
        const double dt_feed = dt * feed;
        const double kill_feed = kill + feed;
        
        for (int y = startY; y < endY; y++) {
            for (int x = 0; x < width; x++) {
                int index = idx(y, x);
                double U = u[index];
                double V = v[index];
                
                double lapU = laplacian(u, y, x);
                double lapV = laplacian(v, y, x);
                
                double reaction = U * V * V;
                double newU = U + dt_Du * lapU - dt * reaction + dt_feed * (1.0 - U);
                double newV = V + dt_Dv * lapV + dt * reaction - dt * kill_feed * V;
                
                u_next[index] = (newU < 0.0) ? 0.0 : (newU > 1.0) ? 1.0 : newU;
                v_next[index] = (newV < 0.0) ? 0.0 : (newV > 1.0) ? 1.0 : newV;
            }
        }
    }
    
    void step() {
        std::vector<std::thread> threads;
        int rowsPerThread = height / numThreads;
        
        for (unsigned int t = 0; t < numThreads; t++) {
            int startY = t * rowsPerThread;
            int endY = (t == numThreads - 1) ? height : (t + 1) * rowsPerThread;
            threads.emplace_back(&TuringPattern::stepRange, this, startY, endY);
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
        
        std::swap(u, u_next);
        std::swap(v, v_next);
    }
    
    void simulate(int steps) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Simulating " << steps << " iterations with " << numThreads << " threads...\n";
        
        for (int i = 0; i < steps; i++) {
            step();
            if ((i + 1) % 1000 == 0) {
                std::cout << "Step " << (i + 1) << "/" << steps << "\r" << std::flush;
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "\nSimulation complete in " << duration.count() / 1000.0 << " seconds!\n";
        std::cout << "Performance: " << (steps * width * height) / (duration.count() / 1000.0) / 1e6 
                  << " million cells/second\n";
    }
    
    void savePPM() {
        // Find min/max for normalization
        double minVal = v[0], maxVal = v[0];
        size_t size = width * height;
        for (size_t i = 0; i < size; i++) {
            minVal = std::min(minVal, v[i]);
            maxVal = std::max(maxVal, v[i]);
        }
        
        // Write PPM
        std::string ppmFile = filename + ".ppm";
        std::ofstream file(ppmFile);
        file << "P3\n" << width << " " << height << "\n255\n";
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int index = idx(y, x);
                int intensity = static_cast<int>(255 * (v[index] - minVal) / (maxVal - minVal));
                int r, g, b;
                if (intensity < 85) {
                    r = 0;
                    g = intensity * 3;
                    b = 255;
                } else if (intensity < 170) {
                    r = (intensity - 85) * 3;
                    g = 255;
                    b = 255 - (intensity - 85) * 3;
                } else {
                    r = 255;
                    g = 255 - (intensity - 170) * 3;
                    b = 0;
                }
                file << r << " " << g << " " << b << " ";
            }
            file << "\n";
        }
        file.close();
        
        // Write metadata
        std::string txtFile = filename + ".txt";
        std::ofstream meta(txtFile);
        meta << std::fixed << std::setprecision(6);
        meta << "Turing Pattern Parameters\n";
        meta << "========================\n";
        meta << "Diffusion U (Du): " << Du << "\n";
        meta << "Diffusion V (Dv): " << Dv << "\n";
        meta << "Feed rate (f): " << feed << "\n";
        meta << "Kill rate (k): " << kill << "\n";
        meta << "Grid size: " << width << "x" << height << "\n";
        meta << "Timestamp: " << std::time(nullptr) << "\n";
        meta.close();
        
        std::cout << "\nFiles saved:\n";
        std::cout << "  Image: " << ppmFile << "\n";
        std::cout << "  Parameters: " << txtFile << "\n";
    }
    
    void printParameters() {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\n=== Turing Pattern Parameters ===\n";
        std::cout << "Du (Diffusion U): " << Du << "\n";
        std::cout << "Dv (Diffusion V): " << Dv << "\n";
        std::cout << "Feed rate: " << feed << "\n";
        std::cout << "Kill rate: " << kill << "\n";
        std::cout << "Grid size: " << width << "x" << height << "\n";
        std::cout << "Output filename: " << filename << "\n";
        std::cout << "================================\n\n";
    }
};

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Generate random parameters for interesting patterns
    std::uniform_real_distribution<> du_dist(0.10, 0.20);
    std::uniform_real_distribution<> dv_dist(0.03, 0.10);
    std::uniform_real_distribution<> feed_dist(0.020, 0.070);
    std::uniform_real_distribution<> kill_dist(0.045, 0.065);
    
    double Du = du_dist(gen);
    double Dv = dv_dist(gen);
    double feed = feed_dist(gen);
    double kill = kill_dist(gen);
    
    int width = 512;
    int height = 512;
    int iterations = 10000;
    
    TuringPattern pattern(width, height, Du, Dv, feed, kill);
    pattern.printParameters();
    
    pattern.simulate(iterations);
    pattern.savePPM();
    
    std::cout << "\nPattern generation complete!\n";
    std::cout << "To convert to PNG: convert <filename>.ppm <filename>.png\n";
    
    return 0;
}
