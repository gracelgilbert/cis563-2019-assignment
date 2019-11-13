#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"

bool selectInGeometry(TV samplePos) {
    return true;
}

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up particle system
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

    // Random sampling of points in a cube or mesh
    int N = 64; // Make it even!!
    int N_points = N*N*N;
    T dx = (T) 1 / ((T)N - (T)1);
    m.resize(N_points);
    x.resize(N_points);
    v.resize(N_points);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                TV samplePos = TV::Zero();
                samplePos(0) = (i-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                samplePos(1) = (i-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                samplePos(2) = (i-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                if (selectInGeometry(samplePos)) {
                    int id = i + j * (N - 1) + k * (N - 1) * (N - 1);
                    m[id] = (T)1/N_points;
                    x[id] = samplePos;
                    v[id] = TV::Zero();
                }

            }
        }
    }

    // set up grid system
    std::vector<T> grid_m;
    std::vector<TV> grid_v;
    std::vector<TV> grid_vn;
    std::vector<TV> grid_f;
    std::vector<int> active_nodes;
    int grid_N = N / 2;
    float grid_h = 1.0 / ((float)grid_N - 1.0);
    // std::cout << grid_N << std::endl;
    int N_cells = grid_N * grid_N * grid_N;
    grid_m.resize(N_cells);
    grid_v.resize(N_cells);
    grid_vn.resize(N_cells);
    grid_f.resize(N_cells);

    // simulate
    driver.ps.m = m;
    driver.ps.v = v;
    driver.ps.x = x;

    driver.gs.res = grid_N;
    driver.gs.h = grid_h;
    driver.gs.m = grid_m;
    driver.gs.v = grid_v;
    driver.gs.vn = grid_vn;
    driver.gs.f = grid_f;
    driver.gs.active_nodes = active_nodes;

    driver.run(24);

    return 0;
}
