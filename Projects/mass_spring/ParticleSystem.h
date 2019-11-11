#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class ParticleSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

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
