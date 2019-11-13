#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class GridSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<T> m;
    std::vector<TV> v;
    std::vector<TV> vn;
    std::vector<TV> f;
    std::vector<int> active_nodes;
    int res;
    float h;

    GridSystem()
    {}
};
