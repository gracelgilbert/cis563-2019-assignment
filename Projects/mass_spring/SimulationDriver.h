#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "ParticleSystem.h"
#include "GridSystem.h"

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
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
      : dt((T)0.0015)  // 150 times bigger dt than explicit. We can't go arbitrarily large because we are still doing approximations to the non-linear problem using taylor expansion.
    {
        gravity = TV::Zero();
        gravity(1) = -9.8;

        sphere_center = TV::Ones()*0.5;
        sphere_radius = 0.2;
        ground = 0.1;
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

        // update velocity
        updateGridVelocity();

        // boundary conditions
        setBoundaryVelocities(3);


        // G2P (including particle advection)
        // Lg = computeGridMomentum();
        // std::cout << "Particle Momentum before g2p: " << Lg[0] << ", " << Lg[1] << ", " << Lg[2] << std::endl;

        // evolveF();
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

    // void evolveF() {

    // }

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



    

    // void advanceOneStepExplicitIntegration()
    // {
    //     int N_points = ms.x.size();
    //     int N_dof = dim*N_points;
	// std::vector<TV> f_spring;
    //     ms.evaluateSpringForces(f_spring);
	// std::vector<TV> f_damping;
	// ms.evaluateDampingForces(f_damping);
	
	// for(int p=0; p<N_points; p++){
    //         if(ms.node_is_fixed[p]){
	//       ms.v[p] = TV::Zero();
    //         }
	//     else{
	//       ms.v[p] += ((f_spring[p]+f_damping[p])/ms.m[p]+gravity)*dt;
	//       ms.x[p] += ms.v[p]*dt;
	//     }
    //     }
    // }
    
    // void advanceOneStepImplicitIntegration()
    // {
    //     int N_points = ms.x.size();
    //     int N_dof = dim*N_points;
    //     SpMat A(N_dof,N_dof);
    //     A.reserve(Eigen::VectorXi::Constant(N_dof,dim*20)); // estimate non-zero entries per column

    //     // build right hand side
    //     std::vector<TV> fn;
    //     ms.evaluateSpringForces(fn);

	// Vec rhs(N_dof);
    //     rhs.setZero();
	// ///////////////////////////////////////////////
	// //  ASSIGNMENT ////////////////////////////////
	// //  Add you code to build f /////////////////// 
	// ///////////////////////////////////////////////

    //     for(int p=0; p<N_points; p++) {
    //         TV vel = ms.v[p];
    //         TV springForce = fn[p];
    //         T mass = ms.m[p];
    //         for (int d = 0; d < dim; d++) {
    //             int i = p * dim + d;
    //             rhs[i] += mass * vel[d] / dt;
    //             rhs[i] += springForce[d];
    //             rhs[i] += mass * gravity[d];
    //         }
    //     }
  
    //     // build the matrix
    //     // Mass matrix contribution (assembly to global)
    //     for(int p=0; p<N_points; p++) {
    //         for (int d = 0; d < dim; d++) {
    //             int i = p * dim + d; // global dof index
    //             A.coeffRef(i, i) += ms.m[p] / (dt * dt);
    //         }
    //     }

    //     for(size_t e=0; e<ms.segments.size(); e++)
    //     {
    //         int particle[2]; // global particle index
    //         particle[0] = ms.segments[e](0);
    //         particle[1] = ms.segments[e](1);


    //         T l0 = ms.rest_length[e];
    //         T E = ms.youngs_modulus;
    //         T l = (ms.x[particle[0]]-ms.x[particle[1]]).norm();
    //         TV n = (ms.x[particle[0]]-ms.x[particle[1]])/l;
    //         T b = ms.damping_coeff;

    //         // Damping matrix contribution
    //         Eigen::Matrix<T,dim,dim> bnn = b * n * n.transpose();
    //         Eigen::Matrix<T,dim*2,dim*2> G_local;
    //         G_local.template block<dim,dim>(0,0) = -bnn;
    //         G_local.template block<dim,dim>(dim,0) = bnn;
    //         G_local.template block<dim,dim>(0,dim) = bnn;
    //         G_local.template block<dim,dim>(dim,dim) = -bnn;

	//     ////////////////////////////////////////////////////////////////// 
	//     //  ASSIGNMENT /////////////////////////////////////////////////// 
	//     //  Add you code to construct local elasticity matrix K_local
	//     /////////////////////////////////////////////// //////////////////
    //         Eigen::Matrix<T,dim,dim> I = Eigen::Matrix<T, dim, dim>::Identity();
    //         Eigen::Matrix<T,dim,dim> K = E * ((1 / l0) - (1 / l)) * (I - n * n.transpose()) + (E / l0) * n * n.transpose();
    //         Eigen::Matrix<T,dim*2,dim*2> K_local;
    //         K_local.template block<dim,dim>(0,0) = -K;
    //         K_local.template block<dim,dim>(dim,0) = K;
    //         K_local.template block<dim,dim>(0,dim) = K;
    //         K_local.template block<dim,dim>(dim,dim) = -K;


	//     ////////////////////////////////////////////////////////////////// 
	//     //  ASSIGNMENT /////////////////////////////////////////////////// 
	//     //  Add you code to add contributions of elasticity and damping to A
	//     //  Note that you need to take care of dirichlet-0 nodes in the
	//     //     corresponding row and columns (by keeping those entries 0)
	//     /////////////////////////////////////////////// //////////////////

    //         for (int i = 0; i < 2; i++) {
    //             for (int j = 0; j < 2; j++) {

    //                 for (int x = 0; x < dim; x++) {
    //                     for (int y = 0; y < dim; y++) {

    //                         if (!ms.node_is_fixed[particle[i]] && !ms.node_is_fixed[particle[j]]) {
    //                             int xindex = dim * particle[i] + x;
    //                             int yindex = dim * particle[j] + y;

    //                             int localxIndex = dim * i + x;
    //                             int localyIndex = dim * j + y;

    //                             A.coeffRef(xindex, yindex) -= G_local(localxIndex, localyIndex) / dt;
    //                             A.coeffRef(xindex, yindex) -= K_local(localxIndex, localyIndex);
    //                         }

    //                     }
    //                 }

    //             }
    //         }
    //     }

    //     // process dirichlet-0 nodes at the diagonal of A
    //     for (size_t p=0; p<ms.node_is_fixed.size(); p++){
    //         if(ms.node_is_fixed[p]){
    //             for(int d=0; d<dim; d++) {
    //                 A.coeffRef(dim * p + d, dim * p + d) = 1;
    //                 rhs[p*dim+d] = 0;
    //             }
    //         }
    //     }

    //     // Eigen::ConjugateGradient<SpMat , Eigen::Lower|Eigen::Upper> krylov;
    //     Eigen::MINRES<SpMat, Eigen::Lower|Eigen::Upper> krylov;

    //     krylov.setTolerance((T)1e-7);
    //     krylov.compute(A);
    //     Vec dx = krylov.solve(rhs);
    //     std::cout << "#iterations:     " << krylov.iterations() << std::endl;
    //     std::cout << "estimated error: " << krylov.error()      << std::endl;

    //     for(int p=0; p<N_points; p++){

    //         TV tentative_dx;
    //         for(int d=0; d<dim; d++)
    //             tentative_dx(d) = dx(dim*p+d);

    //         TV x0 = ms.x[p];
    //         TV v0 = tentative_dx/dt;
    //         TV xt = x0 + tentative_dx;

    //         // ground impulse
    //         if(xt(1) < ground){
    //             T v0n = (ground - x0(1))/dt;

	// 	// This is some hack friction. Real friction force should depend on normal force magnitude.
	// 	T friction_coefficient = 0.1;
	// 	v0(0)*=friction_coefficient;
	// 	v0(2)*=friction_coefficient;
		
    //             v0(1) = v0n;
    //             xt = x0 + v0*dt;
    //             tentative_dx = xt-x0;
    //         }

    //         // sphere impulse
    //         T distance_to_sphere_center = (xt-sphere_center).norm();
    //         if (distance_to_sphere_center < sphere_radius)
    //         {
    //             TV inward_normal = (sphere_center-xt).normalized();
    //             TV v0normal = v0.dot(inward_normal) * inward_normal;
    //             TV v0tangent = v0 - v0normal;
    //             v0normal -= (sphere_radius-distance_to_sphere_center)/dt*inward_normal;
    //             v0 = v0normal + v0tangent;
    //             xt = x0 + v0*dt;
    //             tentative_dx = xt-x0;
    //         }

    //         ms.x[p] += tentative_dx;
    //         ms.v[p] = tentative_dx/dt;
    //     }

    // }



};
