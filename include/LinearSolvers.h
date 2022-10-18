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
       

        //static int constrained_solution(vec_d &A, vec_d &b, std::vector<vec_d> &C, vec_d &d);

    private:
        std::vector<vec_i> fragments;
        std::chrono::steady_clock::duration wtime_fragment, wtime_total;
        std::chrono::steady_clock::duration ptime_fragment, ptime_total;

        void solve_restrained(vec_d &A_mat, vec_d &rhs_in, vec_d &guess, int n_sites, double rho_max, vec_d &out);
        void apply_constraints(vec_d &A_mat, vec_d &rhs, std::vector<vec_d> &constraints, vec_d &constraint_vals);
        void project_constraints(vec_d &metric, vec_d &guess, std::vector<vec_d> &constraints, vec_d &constraint_vals, double max_abs_val);
        inline double sign(double x);
};
#endif // LINEAR_SOLVER