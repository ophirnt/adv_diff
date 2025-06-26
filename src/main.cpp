#include <iostream>
#include <vector>
#include "Mesh.hpp"
#include "Solvers1D.hpp"
#include "IO.hpp"

// Vec assembleField(SimpleMesh *mesh) {
//     Vec T;
//     VecCreateSeq(PETSC_COMM_WORLD, mesh->nx(), &T);

//     for (int i = 0; i < mesh->nx() - 1; i++) {
//         VecSetValue(T, i, 0.0, INSERT_VALUES);
//     }
//     VecSetValue(T, mesh->nx() - 1, 1.0, INSERT_VALUES);
//     VecAssemblyBegin(T);
//     VecAssemblyEnd(T);

//     return T;
// }





int main(int argc, char* argv[]) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    double const    u  = 1.0;
    double const    D  = 0.5;
    double const    x0 = 0;
    double const    x1 = 1.0;
    int    const    nx = 10;
    double const    T0 = 0;
    double const    Tf = 1;

    SimpleMesh mesh = SimpleMesh(0, 1, 9);

    double const Pe = u*mesh.dx() / (2*D);

    std::cout << "Pe = " << Pe << std::endl;

    std::vector<double> T(mesh.nx());
    T[0] = T0;
    T[mesh.nx() - 1] = Tf;
    print_vector(T);

    solve_upwind(T, &mesh, 1.0, 0.5);

    std::cout << "Solved with upwind" << std::endl;
    print_vector(T);

    solve_cds(T, &mesh, 1.0, 0.5);
    print_vector(T);

    std::cout << "Exact solution" << std::endl;
    std::vector<double> Texact(mesh.nx());
    Texact[mesh.nx() - 1] = 1.0;
    solve_exact(Texact, &mesh, u, D);
    print_vector(Texact);
    mesh.print_mesh();


    PetscFinalize();
    return 0;
}





