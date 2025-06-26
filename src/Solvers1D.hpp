#include <iostream>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <vector>
#include "Mesh.hpp"


#ifndef SOLVERS1D_HPP
#define SOLVERS1D_HPP

void solve_cds(std::vector<double> &T, SimpleMesh *mesh, double u, double D) {

    const double adv = u * 1.0 / 2;
    const double diff = D / mesh->dx() * 1.0;
    const double west_adv = (std::abs(adv) + adv) / 2.0;
    const double east_adv = (std::abs(adv) - adv) / 2.0;

    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, mesh->nx() - 2, mesh->nx() - 2, mesh->nx() - 2, mesh->nx() - 2);
    MatSetType(A, MATSEQAIJ);
    MatSetUp(A);

    for (int i = 0; i < mesh->nx() - 2; ++i) {
        if (i > 0)
            MatSetValue(A, i, i - 1, -diff -adv, INSERT_VALUES);

        MatSetValue(A, i, i, 2 * diff, INSERT_VALUES);

        if (i < mesh->nx() - 3)
            MatSetValue(A, i, i + 1, -diff +adv, INSERT_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatView(A, PETSC_VIEWER_STDOUT_SELF);

    Vec b;
    MatCreateVecs(A, &b, NULL);

    VecSetValue(b, 0, T[0] * (diff - west_adv), INSERT_VALUES);
    VecSetValue(b, mesh->nx() - 3, T[mesh->nx() - 1] * (diff - east_adv), INSERT_VALUES);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    VecView(b, PETSC_VIEWER_STDOUT_SELF);

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);  // matrix and preconditioner
    KSPSetFromOptions(ksp);     // allow runtime flags

    Vec T_internal;
    VecCreate(PETSC_COMM_SELF, &T_internal);
    VecSetSizes(T_internal, mesh->nx()-2, mesh->nx()-2);
    VecSetFromOptions(T_internal);

    // Solve Ax = b
    KSPSolve(ksp, b, T_internal);

    // Optional: View solution and solver stats
    KSPView(ksp, PETSC_VIEWER_STDOUT_SELF);
    VecView(T_internal, PETSC_VIEWER_STDOUT_SELF);

    PetscScalar temp;
    for(int i = 1; i < mesh->nx() - 1; i++){
        PetscInt idx1 = i - 1;  // T_internal index
        VecGetValues(T_internal, 1, &idx1, &temp);
        T[i] = temp;
        //VecSetValue(T, i, temp, INSERT_VALUES);
    }

    // Cleanup
    VecDestroy(&b);
    KSPDestroy(&ksp);
    VecDestroy(&T_internal);
}


void solve_upwind(std::vector<double> &T, SimpleMesh *mesh, double u, double D) {
    
    const double adv = u * 1.0;
    const double diff = D / mesh->dx() * 1.0;
    const double west_adv = (std::abs(adv) + adv) / 2.0;
    const double east_adv = (std::abs(adv) - adv) / 2.0;

    std::cout << "PEKO" << adv << "," << diff << "," << west_adv << "," << east_adv << std::endl;

    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, mesh->nx() - 2, mesh->nx() - 2, mesh->nx() - 2, mesh->nx() - 2);
    MatSetType(A, MATSEQAIJ);
    MatSetUp(A);

    for (int i = 0; i < mesh->nx() - 2; ++i) {
        if (i > 0)
            MatSetValue(A, i, i - 1, -diff + west_adv, INSERT_VALUES);

        MatSetValue(A, i, i, 2 * diff - adv, INSERT_VALUES);

        if (i < mesh->nx() - 3)
            MatSetValue(A, i, i + 1, -diff + east_adv, INSERT_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatView(A, PETSC_VIEWER_STDOUT_SELF);

    Vec b;
    MatCreateVecs(A, &b, NULL);

    VecSetValue(b, 0, T[0] * (diff - west_adv), INSERT_VALUES);
    VecSetValue(b, mesh->nx() - 3, T[mesh->nx() - 1] * (diff - east_adv), INSERT_VALUES);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);  // matrix and preconditioner
    KSPSetFromOptions(ksp);     // allow runtime flags

    Vec T_internal;
    VecCreate(PETSC_COMM_SELF, &T_internal);
    VecSetSizes(T_internal, mesh->nx()-2, mesh->nx()-2);
    VecSetFromOptions(T_internal);

    // Solve Ax = b
    KSPSolve(ksp, b, T_internal);

    // Optional: View solution and solver stats
    KSPView(ksp, PETSC_VIEWER_STDOUT_SELF);
    VecView(T_internal, PETSC_VIEWER_STDOUT_SELF);

    PetscScalar temp;
    for(int i = 1; i < mesh->nx() - 1; i++){
        PetscInt idx1 = i - 1;  // T_internal index
        VecGetValues(T_internal, 1, &idx1, &temp);
        T[i] = temp;
        //VecSetValue(T, i, temp, INSERT_VALUES);
    }

    // Cleanup
    VecDestroy(&b);
    KSPDestroy(&ksp);
    VecDestroy(&T_internal);
}

void solve_exact(std::vector<double> &T, SimpleMesh *mesh, double u, double D){
    double const k = exp(u/D * mesh->x(mesh->nx() - 1));
    double const C2 = (T[0] - T[mesh->nx() - 1]) / (1 - k);
    double const C1 = T[0] - C2;


    for(int i = 0; i < mesh->nx(); i++){
        double solution = C1 + C2 * exp(u/D * mesh->x(i));
        T[i] = solution;
    }
}

#endif