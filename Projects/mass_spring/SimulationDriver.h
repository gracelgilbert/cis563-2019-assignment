#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "ParticleSystem.h"
#include "GridSystem.h"
#include <Eigen/SVD>

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    ParticleSystem<T,dim> ps;
    GridSystem<T, dim> gs;
    T dt;
    TV gravity;

    TV sphere_center;
    T sphere_radius;
    T ground;

    SimulationDriver()
      // : dt((T)0.00001) 
      : dt((T)0.0005)  // 150 times bigger dt than explicit. We can't go arbitrarily large because we are still doing approximations to the non-linear problem using taylor expansion.
    {
        gravity = TV::Zero();
        gravity(1) = -9.8;

        sphere_center = TV::Ones()*0.5;
        sphere_radius = 0.2;
        ground = 0.1;
    }

    bool selectInGeometry(TV samplePos) {

        return (samplePos[0] > 0.3 && samplePos[0] < 0.7 && 
                samplePos[1] > 0.5 && samplePos[1] < 0.9 &&
                samplePos[2] > 0.3 && samplePos[2] < 0.7);
    }

    void run(const int max_frame)
    {
        mkdir("output/", 0777);
        std::string filename = "output/" + std::to_string(max_frame) + ".poly";
        ps.dumpPoly(filename);
        std::cout << std::endl;

        
        for(int frame=1; frame<=max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                std::cout << "Step " << step << std::endl;
                advanceOneStepGridUpdate();
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            ps.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceOneStepGridUpdate()
    {
        // Zero out grid
        gs.m.clear();
        gs.v.clear();
        gs.vn.clear();
        gs.f.clear(); // TODO: are the clears necessary?
        gs.active_nodes.clear();
        gs.m.resize(gs.res * gs.res * gs.res, (T)0);
        gs.v.resize(gs.res * gs.res * gs.res, TV::Zero());
        gs.vn.resize(gs.res * gs.res * gs.res, TV::Zero());
        gs.f.resize(gs.res * gs.res * gs.res, TV::Zero());

        // TV Lp = computeParticleMomentum();
        // std::cout << "Particle Momentum before p2g: " << Lp[0] << ", " << Lp[1] << ", " << Lp[2] << std::endl;

        transferP2G();


        // TV Lg = computeGridMomentum();
        // std::cout << "Particle Momentum after p2g: " << Lg[0] << ", " << Lg[1] << ", " << Lg[2] << std::endl;

        // compute force
        addGravity();
        addElasticity();

        // update velocity
        updateGridVelocity();

        // boundary conditions
        setBoundaryVelocities(2);



        // G2P (including particle advection)
        // Lg = computeGridMomentum();
        // std::cout << "Particle Momentum before g2p: " << Lg[0] << ", " << Lg[1] << ", " << Lg[2] << std::endl;

        evolveF();
        transferG2P(0.95);
    }

    TV computeParticleMomentum() {

        TV totalMomentum = TV::Zero();
        for (int p = 0; p < (int)ps.m.size(); ++p) {
            totalMomentum += ps.m[p] * ps.v[p];
        }
        return totalMomentum;
    }

    void transferP2G() {
        for (int p = 0; p < (int)ps.m.size(); ++p) {
            TV X = ps.x[p];
            TV X_index_space = X / gs.h;
            int base_node1 = computeBaseNode(X_index_space[0]);
            int base_node2 = computeBaseNode(X_index_space[1]);
            int base_node3 = computeBaseNode(X_index_space[2]);
            
            TV w1 = ComputeWeights1D(X_index_space[0], base_node1);
            TV w2 = ComputeWeights1D(X_index_space[1], base_node2);
            TV w3 = ComputeWeights1D(X_index_space[2], base_node3);
   

            for (int i = 0; i < dim; ++i) {
                float wi = w1[i];
                int node_i = base_node1 + (i);
                for (int j = 0; j < dim; ++j) {
                    float wij = wi * w2[j];
                    int node_j = base_node2 + (j);
                    for (int k = 0; k < dim; ++k) {
                        float wijk = wij * w3[k];
                        int node_k = base_node3 + (k);

                        // splat mass
                        int gridIndex = node_i + node_j * gs.res + node_k * gs.res * gs.res;
                        gs.m[gridIndex] = gs.m[gridIndex] + ps.m[p] * wijk;
                        // splat momentum
                        for (int d = 0; d < dim; ++d) {
                            gs.vn[gridIndex][d] = gs.vn[gridIndex][d] + (wijk * ps.m[p]) * ps.v[p][d];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < (int)gs.m.size(); ++i) {
            if (gs.m[i] != 0) {
                gs.active_nodes.push_back(i);
                gs.vn[i] = gs.vn[i] /  gs.m[i];
            } else {
                gs.vn[i] = TV::Zero();
            }
        }
    }

    int computeBaseNode(float x) {
        return floor(x - 0.5) + 1;
    }

    TV ComputeWeights1D(float x, int base_node) {
        TV w = TV::Zero();

        float d0 = x - base_node + 1;
        float z = 1.5 - d0;
        float z2 = z * z;
        w[0] = 0.5 * z2;

        float d1 = d0 - 1;
        w[1] = 0.75 - d1 * d1;

        float d2 = 1 - d1;
        float zz = 1.5 - d2;
        float zz2 = zz * zz;
        w[2] = 0.5 * zz2;
        return w;
    }

    void ComputeWeightsWithGradients1D(float x, int base_node, TV& w, TV& dw) {
        w = TV::Zero();
        dw = TV::Zero();

        float d0 = x - base_node + 1;
        float z = 1.5 - d0;
        float z2 = z * z;
        w[0] = 0.5 * z2;

        float d1 = d0 - 1;
        w[1] = 0.75 - d1 * d1;

        float d2 = 1 - d1;
        float zz = 1.5 - d2;
        float zz2 = zz * zz;
        w[2] = 0.5 * zz2;

        dw[0] = -z;
        dw[1] = -2 * d1;
        dw[2] = zz;
    }

    TV computeGridMomentum() {
        TV totalMomentum = TV::Zero();
        for (int i = 0; i < (int)gs.m.size(); ++i) {
            totalMomentum += gs.m[i] * gs.vn[i];
        }
    }

    void addGravity() {
        for (int i = 0; i < (int)gs.active_nodes.size(); ++i) {
            int currIndex = gs.active_nodes[i];
            gs.f[currIndex] = gs.f[currIndex] + gs.m[currIndex] * gravity;
        }
    }

    void polarSVD(TM F, TM& U, TM& Sigma, TM& V) {
        Eigen::JacobiSVD<TM> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
        U = svd.matrixU();
        V = svd.matrixV();
        Sigma = svd.singularValues().asDiagonal();
        float detU = U.determinant();
        float detV = V.determinant();
        if (detU < 0) {
            U(0, 2) = U(0, 2) * -1;
            U(1, 2) = U(1, 2) * -1;
            U(2, 2) = U(2, 2) * -1;
            Sigma(2, 2) = Sigma(2, 2) * -1;
        }
        if (detV < 0) {
            V(0, 2) = V(0, 2) * -1;
            V(1, 2) = V(1, 2) * -1;
            V(2, 2) = V(2, 2) * -1;
            Sigma(2, 2) = Sigma(2, 2) * -1;
        }
    }

    TM fixedCorotated(TM F) {
        TM U;
        TM Sigma;
        TM V;
        polarSVD(F, U, Sigma, V);
        
        TM R = U * V.transpose();
        TM tempF = U * Sigma * V.transpose();
        float J = tempF.determinant();

        TM A = (tempF.adjoint()).transpose();
        // A(0, 0) =   F(2, 2) * F(1, 1) - F(2, 1) * F(1, 2);
        // A(1, 0) = -(F(2, 2) * F(0, 1) - F(2, 1) * F(0, 2));
        // A(2, 0) =   F(1, 2) * F(0, 1) - F(1, 1) * F(0, 2);
        // A(0, 1) = -(F(2, 2) * F(1, 0) - F(2, 0) * F(1, 2));
        // A(1, 1) =   F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
        // A(2, 1) = -(F(1, 2) * F(0, 0) - F(1, 0) * F(0, 2));
        // A(0, 2) =   F(2, 1) * F(1, 0) - F(2, 0) * F(1, 1);
        // A(1, 2) = -(F(2, 1) * F(0, 0) - F(2, 0) * F(0, 1));
        // A(2, 2) =   F(1, 1) * F(0, 0) - F(1, 0) * F(0, 1);

        return 2 * ps.mu * (tempF - R) + ps.lambda * (J - 1) * A;
    }

    void addElasticity() {
        for (int p = 0; p < (int)ps.m.size(); ++p) {
            TM F = ps.F[p];
            TM P = fixedCorotated(F);
            TM V0PFt = ps.V0[p] * P * F.transpose();
            // std::cout << ps.V0[p] << std::endl;

            TV X = ps.x[p];
            TV X_index_space = X / gs.h;
            int base_node1 = computeBaseNode(X_index_space[0]);
            int base_node2 = computeBaseNode(X_index_space[1]);
            int base_node3 = computeBaseNode(X_index_space[2]);

            TV w1;
            TV w2;
            TV w3;
            TV dw1;
            TV dw2;
            TV dw3;

            ComputeWeightsWithGradients1D(X_index_space[0], base_node1, w1, dw1);
            ComputeWeightsWithGradients1D(X_index_space[1], base_node2, w2, dw2);
            ComputeWeightsWithGradients1D(X_index_space[2], base_node3, w3, dw3);

            for (int i = 0; i < dim; ++i) {
                float wi = w1[i];
                float dwidxi = dw1[i] / gs.h;

                int node_i = base_node1 + (i);

                for (int j = 0; j < dim; ++j) {
                    float wj = w2[j];
                    float wij = wi * wj;
                    float dwijdxi = dwidxi * wj;
                    float dwijdxj = wi / gs.h * dw2[j];

                    int node_j = base_node2 + (j);

                    for (int k = 0; k < dim; ++k) {
                        float wijk = wij * w3[k];
                        float dwijkdxi = dwijdxi * w3[k];
                        float dwijkdxj = dwijdxj * w3[k];
                        float dwijkdxk = wij / gs.h * dw3[k];

                        int node_k = base_node3 + (k);

                        TV grad_w = TV::Zero();
                        grad_w[0] = dwijkdxi;
                        grad_w[1] = dwijkdxj;
                        grad_w[2] = dwijkdxk;

                        TV addForce = -V0PFt * grad_w;

                        int currIndex = node_i + node_j * gs.res + node_k * gs.res * gs.res;
                        gs.f[currIndex] = gs.f[currIndex] + addForce;
                    }
                }
            }

        }
    }

    void evolveF() {
        for (int p = 0; p < (int)ps.m.size(); ++p) {
            TM thisF = ps.F[p];

            TV X = ps.x[p];
            TV X_index_space = X / gs.h;

            int base_node1 = computeBaseNode(X_index_space[0]);
            int base_node2 = computeBaseNode(X_index_space[1]);
            int base_node3 = computeBaseNode(X_index_space[2]);

            TV w1;
            TV w2;
            TV w3;
            TV dw1;
            TV dw2;
            TV dw3;

            ComputeWeightsWithGradients1D(X_index_space[0], base_node1, w1, dw1);
            ComputeWeightsWithGradients1D(X_index_space[1], base_node2, w2, dw2);
            ComputeWeightsWithGradients1D(X_index_space[2], base_node3, w3, dw3);

            TM grad_vp = TM::Zero();
            for (int i = 0; i < dim; ++i) {
                float wi = w1[i];
                float dwidxi = dw1[i] / gs.h;

                int node_i = base_node1 + (i);

                for (int j = 0; j < dim; ++j) {
                    float wj = w2[j];
                    float wij = wi * wj;
                    float dwijdxi = dwidxi * wj;
                    float dwijdxj = wi / gs.h * dw2[j];

                    int node_j = base_node2 + (j);

                    for (int k = 0; k < dim; ++k) {
                        float wijk = wij * w3[k];
                        float dwijkdxi = dwijdxi * w3[k];
                        float dwijkdxj = dwijdxj * w3[k];
                        float dwijkdxk = wij * dw3[k] / gs.h;

                        int node_k = base_node3 + (k);

                        TV grad_w = TV::Zero();
                        grad_w[0] = dwijkdxi;
                        grad_w[1] = dwijkdxj;
                        grad_w[2] = dwijkdxk;

                        int currIndex = node_i + node_j * gs.res + node_k * gs.res * gs.res;
                        if (currIndex >= (int) gs.v.size() - 1) {
                            std::cout << "ERROR: Grid out of bounds" << std::endl;
                        }

                        TV vijk = gs.v[currIndex];
                        grad_vp = grad_vp + vijk * grad_w.transpose();
                    }
                }
            }

            TM newFP = (TM::Identity(dim, dim) + dt * grad_vp) * thisF;
            // std::cout << newFP << std::endl;
            // std::cout << std::endl;
            // for (int i = 0; i < dim; ++i) {
            //     for (int j = 0; j < dim; ++j) {
            //         ps.F[p](i, j) = newFP(i, j);
            //     }
            // }
            ps.F[p] = newFP;
        }
    }


    void updateGridVelocity(){
        for (int i = 0; i < (int)gs.active_nodes.size(); ++i) {
            int currIndex = gs.active_nodes[i];
            gs.v[currIndex] = gs.vn[currIndex] + dt * gs.f[currIndex] / gs.m[currIndex];
        }
    }

    void setBoundaryVelocities(float thickness) {
        // int Nx = gs.v.size();
        // int Ny = dim;
        // int Nz = dim;

        for (int i = 0; i < thickness; ++i) {
            for (int j = 0; j < gs.res; ++j) {
                for (int k = 0; k < gs.res; ++k) {
                    int gridIndex = i + j * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        for (int i = gs.res - thickness + 1; i < gs.res; ++i) {
            for (int j = 0; j < gs.res; ++j) {
                for (int k = 0; k < gs.res; ++k) {
                    int gridIndex = j + i * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        for (int i = 0; i < gs.res; ++i) {
            for (int j = 0; j < thickness; ++j) {
                for (int k = 0; k < gs.res; ++k) {
                    int gridIndex = i + j * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        for (int i = 0; i < gs.res; ++i) {
            for (int j = gs.res - thickness + 1; j < gs.res; ++j) {
                for (int k = 0; k < gs.res; ++k) {
                    int gridIndex = j + i * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        for (int i = 0; i < gs.res; ++i) {
            for (int j = 0; j < gs.res; ++j) {
                for (int k = 0; k < thickness; ++k) {
                    int gridIndex = i + j * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        for (int i = 0; i < gs.res; ++i) {
            for (int j = 0; j < gs.res; ++j) {
                for (int k = gs.res - thickness; k < gs.res; ++k) {
                    int gridIndex = j + i * gs.res + k * gs.res * gs.res;
                    gs.v[gridIndex] = TV::Zero();
                }
            }
        }

        // for (int j = 0; j < thickness; ++j) {
        //     for (int i = 0; i < gs.res; ++i) {
        //         for (int k = 0; k < gs.res; ++k) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }

        // for (int j = gs.res - thickness - 1; j < gs.res; ++j) {
        //     for (int i = 0; j < gs.res; ++i) {
        //         for (int k = 0; k < gs.res; ++k) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }

        // for (int k = 0; k < thickness; ++k) {
        //     for (int j = 0; j < gs.res; ++j) {
        //         for (int i = 0; i < gs.res; ++i) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }

        // for (int k = gs.res - thickness - 1; k < gs.res; ++k) {
        //     for (int j = 0; j < gs.res; ++j) {
        //         for (int i = 0; i < gs.res; ++i) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }


        // for (int i = Nx - thickness; i < Nx; ++i) {
        //     for (int j = 0; j < Ny; ++j) {
        //         for (int k = 0; k < Nz; ++k) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }
        // for (int i = 0; i < Nx; ++i) {
        //     for (int j = 0; j < dim; ++j) {
        //         for (int k = 0; k < dim; ++k) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }

        // for (int i = 0; i < Nx; ++i) {
        //     for (int j = Ny - thickness; j < Ny; ++j) {
        //         for (int k = 0; k < Nz; ++k) {
        //             int gridIndex = i + j * gs.res + k * gs.res * gs.res;
        //             gs.v[gridIndex] = TV::Zero();
        //         }
        //     }
        // }


    }

    void transferG2P(float flip) {
        for (int p = 0; p < (int)ps.m.size(); ++p) {
            TV X = ps.x[p];
            TV X_index_space = X / gs.h;
            int base_node1 = computeBaseNode(X_index_space[0]);
            int base_node2 = computeBaseNode(X_index_space[1]);
            int base_node3 = computeBaseNode(X_index_space[2]);
            
            TV w1 = ComputeWeights1D(X_index_space[0], base_node1);
            TV w2 = ComputeWeights1D(X_index_space[1], base_node2);
            TV w3 = ComputeWeights1D(X_index_space[2], base_node3);


            TV vpic = TV::Zero();
            TV vflip = ps.v[p];


            for (int i = 0; i < dim; ++i) {
                float wi = w1[i];
                int node_i = base_node1 + (i);
                for (int j = 0; j < dim; ++j) {
                    float wij = wi * w2[j];
                    int node_j = base_node2 + (j);
                    for (int k = 0; k < dim; ++k) {
                        float wijk = wij * w3[k];
                        int node_k = base_node3 + (k);
                        int gridIndex = node_i + node_j * gs.res + node_k * gs.res * gs.res;

                        vpic = vpic + wijk * gs.v[gridIndex];
                        vflip = vflip + wijk * (gs.v[gridIndex] - gs.vn[gridIndex]);
                    }
                }
            }

            ps.v[p] = (1 - flip) * vpic + flip * vflip;
            ps.x[p] = ps.x[p] + dt * vpic;

        }
    }




};
