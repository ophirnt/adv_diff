#ifndef SIMPLEMESH_HPP
#define SIMPLEMESH_HPP

#include <vector>

class SimpleMesh {
public:
    SimpleMesh(double x0, double x1, int nx);
    ~SimpleMesh() = default;

    double x(int i) const;
    double dx() const;
    int nx() const;
private:
    double dx_;
    std::vector<double> x_;
    int nx_;
};

// Constructor
SimpleMesh::SimpleMesh(double x0, double x1, int nx) : x_(nx + 2), nx_(nx) {
    dx_ = (x1 - x0) / nx;
    x_[0] = x0;
    for (int i = 1; i <= nx; i++) {
        x_[i] = x_[i - 1] + dx_;
    }
}
// Accessor for spacing
double SimpleMesh::dx() const {
    return dx_;
}

// Accessor for position
double SimpleMesh::x(int i) const {
    return x_[i];
}

int SimpleMesh::nx() const {
    return nx_;
}

#endif
