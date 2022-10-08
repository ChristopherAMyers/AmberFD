#include <stdio.h>
#include <vector>
#include <lapacke.h>
#include <cblas.h>
#include <iostream>
#include <cmath>
#include "common.h"
#include <omp.h>
#include <chrono>

#ifndef LINEAR_SOLVER
#define LINEAR_SOLVER

class DivideAndConquer
{
    public:
        DivideAndConquer();
        ~DivideAndConquer();

        void assign_fragments(std::vector<vec_i> &fragments);
        void solve(vec_d &Coulomb_mat, vec_d &pot_vec, std::vector<vec_i> &fragments_in, vec_d &delta_rho_out, vec_d &exponents);

    private:
        std::vector<vec_i> fragments;
        // double wtime_fragment, wtime_total;
        // double ptime_fragment, ptime_total;

        std::chrono::steady_clock::duration wtime_fragment, wtime_total;
        std::chrono::steady_clock::duration ptime_fragment, ptime_total;

};
#endif // LINEAR_SOLVER