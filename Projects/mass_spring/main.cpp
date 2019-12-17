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
#include "mesh_query0.1/mesh_query.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;

    SimulationDriver<T,dim> driver;

    // set up particle system
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<TM> F;
    std::vector<T> V0;


    // Generate Mesh:
    std::vector<tinyobj::shape_t> shapes1;
    std::vector<tinyobj::material_t> materials1;
    tinyobj::attrib_t attribs1;
    std::string warning1;
    std::string error1;
    std::string filename1 = "../obj/alpaca1.obj";
    bool objLoad1 = tinyobj::LoadObj(&attribs1, &shapes1, &materials1, &warning1, &error1, filename1.c_str());
    if (!warning1.empty()) {
        std::cout << "WARN: " << warning1 << std::endl;
    }
    if (!error1.empty()) {
		std::cout << "ERR: " << error1 << std::endl;
	}
    if (!objLoad1) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions1 = std::vector<double>();
    std::vector<int> triangles1 = std::vector<int>();

    for (int i = 0; i < (int) attribs1.vertices.size(); ++i) {
        positions1.push_back((double) attribs1.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes1) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles1.push_back(currIndices.at(i).vertex_index);
        }
    }

    std::vector<tinyobj::shape_t> shapes2;
    std::vector<tinyobj::material_t> materials2;
    tinyobj::attrib_t attribs2;
    std::string warning2;
    std::string error2;
    std::string filename2 = "../obj/alpaca2.obj";
    bool objLoad2 = tinyobj::LoadObj(&attribs2, &shapes2, &materials2, &warning2, &error2, filename2.c_str());
    if (!warning2.empty()) {
        std::cout << "WARN: " << warning2 << std::endl;
    }
    if (!error2.empty()) {
		std::cout << "ERR: " << error2 << std::endl;
	}
    if (!objLoad2) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions2 = std::vector<double>();
    std::vector<int> triangles2 = std::vector<int>();

    for (int i = 0; i < (int) attribs2.vertices.size(); ++i) {
        positions2.push_back((double) attribs2.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes2) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles2.push_back(currIndices.at(i).vertex_index);
        }
    }

    std::vector<tinyobj::shape_t> shapes3;
    std::vector<tinyobj::material_t> materials3;
    tinyobj::attrib_t attribs3;
    std::string warning3;
    std::string error3;
    std::string filename3 = "../obj/alpaca3.obj";
    bool objLoad3 = tinyobj::LoadObj(&attribs3, &shapes3, &materials3, &warning3, &error3, filename3.c_str());
    if (!warning3.empty()) {
        std::cout << "WARN: " << warning3 << std::endl;
    }
    if (!error3.empty()) {
		std::cout << "ERR: " << error3 << std::endl;
	}
    if (!objLoad3) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions3 = std::vector<double>();
    std::vector<int> triangles3 = std::vector<int>();

    for (int i = 0; i < (int) attribs3.vertices.size(); ++i) {
        positions3.push_back((double) attribs3.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes3) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles3.push_back(currIndices.at(i).vertex_index);
        }
    }










    std::vector<tinyobj::shape_t> shapes4;
    std::vector<tinyobj::material_t> materials4;
    tinyobj::attrib_t attribs4;
    std::string warning4;
    std::string error4;
    std::string filename4 = "../obj/torus1.obj";
    bool objLoad4 = tinyobj::LoadObj(&attribs4, &shapes4, &materials4, &warning4, &error4, filename4.c_str());
    if (!warning4.empty()) {
        std::cout << "WARN: " << warning4 << std::endl;
    }
    if (!error4.empty()) {
		std::cout << "ERR: " << error4 << std::endl;
	}
    if (!objLoad4) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions4 = std::vector<double>();
    std::vector<int> triangles4 = std::vector<int>();

    for (int i = 0; i < (int) attribs4.vertices.size(); ++i) {
        positions4.push_back((double) attribs4.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes4) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles4.push_back(currIndices.at(i).vertex_index);
        }
    }

    std::vector<tinyobj::shape_t> shapes5;
    std::vector<tinyobj::material_t> materials5;
    tinyobj::attrib_t attribs5;
    std::string warning5;
    std::string error5;
    std::string filename5 = "../obj/torus2.obj";
    bool objLoad5 = tinyobj::LoadObj(&attribs5, &shapes5, &materials5, &warning5, &error5, filename5.c_str());
    if (!warning5.empty()) {
        std::cout << "WARN: " << warning5 << std::endl;
    }
    if (!error5.empty()) {
		std::cout << "ERR: " << error5 << std::endl;
	}
    if (!objLoad5) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions5 = std::vector<double>();
    std::vector<int> triangles5 = std::vector<int>();

    for (int i = 0; i < (int) attribs5.vertices.size(); ++i) {
        positions5.push_back((double) attribs5.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes5) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles5.push_back(currIndices.at(i).vertex_index);
        }
    }

    std::vector<tinyobj::shape_t> shapes6;
    std::vector<tinyobj::material_t> materials6;
    tinyobj::attrib_t attribs6;
    std::string warning6;
    std::string error6;
    std::string filename6 = "../obj/torus3.obj";
    bool objLoad6 = tinyobj::LoadObj(&attribs6, &shapes6, &materials6, &warning6, &error6, filename6.c_str());
    if (!warning6.empty()) {
        std::cout << "WARN: " << warning6 << std::endl;
    }
    if (!error6.empty()) {
		std::cout << "ERR: " << error6 << std::endl;
	}
    if (!objLoad6) {
		std::cout << "failed to load obj";
	}

    std::vector<double> positions6 = std::vector<double>();
    std::vector<int> triangles6 = std::vector<int>();

    for (int i = 0; i < (int) attribs6.vertices.size(); ++i) {
        positions6.push_back((double) attribs6.vertices.at(i));
    }

    for (tinyobj::shape_t shape : shapes6) {
		// Get the indices of the points in each shape
		std::vector<tinyobj::index_t> &currIndices = shape.mesh.indices;
		// Make sure number of indices is a multiple of 3 for triangulation
		if (currIndices.size() % 3 != 0) {
			std::cout << "not triangles" << std::endl;
		}
        for (int i = 0; i < (int) currIndices.size(); ++i) {
            triangles6.push_back(currIndices.at(i).vertex_index);
        }
    }

    MeshObject* mesh1 = construct_mesh_object(positions1.size(), positions1.data(), (int) (triangles1.size() / 3.0), triangles1.data());
    MeshObject* mesh2 = construct_mesh_object(positions2.size(), positions2.data(), (int) (triangles2.size() / 3.0), triangles2.data());
    MeshObject* mesh3 = construct_mesh_object(positions3.size(), positions3.data(), (int) (triangles3.size() / 3.0), triangles3.data());
    MeshObject* mesh4 = construct_mesh_object(positions4.size(), positions4.data(), (int) (triangles4.size() / 3.0), triangles4.data());
    MeshObject* mesh5 = construct_mesh_object(positions5.size(), positions5.data(), (int) (triangles5.size() / 3.0), triangles5.data());
    MeshObject* mesh6 = construct_mesh_object(positions6.size(), positions6.data(), (int) (triangles6.size() / 3.0), triangles6.data());

    // Random sampling of points in a cube or mesh
    int N = 512; // Make it even!!
    T dx = (T) 1 / ((T)N - (T)1);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                TV samplePos = TV::Zero();
                samplePos(0) = (i-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                samplePos(1) = (j-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                samplePos(2) = (k-1 + ((double) rand() / (RAND_MAX)) + 1)*dx;
                const double point[3] = {samplePos(0), samplePos(1), samplePos(2)};
                if (point_inside_mesh(point, mesh1) || point_inside_mesh(point, mesh2) || point_inside_mesh(point, mesh3) ||
                    point_inside_mesh(point, mesh4) || point_inside_mesh(point, mesh5) || point_inside_mesh(point, mesh6)) {
                    x.push_back(samplePos);
                }
            }
        }
    }
    std::cout << "number of particles: " << x.size() << std::endl;

    float rho = 1000;
    float initialVolume = dx * dx * dx;
    int numParticles = x.size();
    m.resize(numParticles, rho * initialVolume);
    v.resize(numParticles, TV::Zero());
    F.resize(numParticles, TM::Identity(dim, dim));

    // set up grid system
    int grid_N = N / 2;
    float grid_h = 1.0 / ((float)grid_N - 1.0);

    V0.resize(numParticles, grid_h * grid_h * grid_h / 8.0);


    // simulate
    driver.ps.m = m;
    driver.ps.v = v;
    driver.ps.x = x;
    driver.ps.F = F;
    driver.ps.V0 = V0;

    driver.gs.res = grid_N;
    driver.gs.h = grid_h;

    driver.run(240);

    return 0;
}
