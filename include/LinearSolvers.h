#include <stdio.h>
#include <vector>
#include <lapacke.h>
#include <cblas.h>
#include <iostream>
#include <cmath>
#include "common.h"
#include <omp.h>

#ifndef LINEAR_SOLVER
#define LINEAR_SOLVER

class DivideAndConquer
{
    public:
        DivideAndConquer();

        void assign_fragments(std::vector<vec_i> &fragments);
        void solve(vec_d &Coulomb_mat, vec_d &pot_vec, std::vector<vec_i> &fragments_in, vec_d &delta_rho_out, vec_d &exponents);

    private:
        std::vector<vec_i> fragments;

};
#endif // LINEAR_SOLVER