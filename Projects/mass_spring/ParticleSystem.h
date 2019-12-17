#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class ParticleSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;
    
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<TM> F;
    std::vector<T> V0;

    float E = 1e4;
    float nu = 0.2;
    float mu = E / (2 * (1 + nu));
    float lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    float rho = 1000;

    ParticleSystem()
    {}

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        fs << "END\n";
        fs.close();
    }
    
};
