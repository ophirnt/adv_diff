#ifndef SIMPLEMESH_HPP
#define SIMPLEMESH_HPP

#include <vector>

class SimpleMesh {
public:
    SimpleMesh(double x0, double x1, int nx);
    ~SimpleMesh() = default;

    void print_mesh();
    double x(int i) const;
    const std::vector<double>& x() const;
    double dx() const;
    int nx() const;
private:
    double dx_;
    std::vector<double> x_;
    int nx_;
};

// Constructor
SimpleMesh::SimpleMesh(double x0, double x1, int nx) : x_(nx + 2), nx_(nx + 2) {
    dx_ = (x1 - x0) / nx;

    std::cout << dx_ << nx_ << nx << std::endl;

    x_[0] = x0;
    x_[1] = x0 + dx_ / 2;
    for (int i = 2; i < nx_; i++) {
        x_[i] = x_[i - 1] + dx_;
    }
    x_[nx_ - 1] = x_[nx_ - 2] + dx_/2;
}

void SimpleMesh::print_mesh(){
    std::cout << "MESH: ";
    std::cout << x_[0];
    for (int i = 1; i < nx_; i++){
        std::cout << ", " << x_[i]; 
    }
    std::cout << std::endl;
}

// Accessor for spacing
double SimpleMesh::dx() const {
    return dx_;
}

// Accessor for position
double SimpleMesh::x(int i) const {
    return x_[i];
}

const std::vector<double>& SimpleMesh::x() const{
    return x_;
}

int SimpleMesh::nx() const {
    return nx_;
}

#endif
