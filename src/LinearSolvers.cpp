#include "LinearSolvers.h"
using namespace std::chrono;
DivideAndConquer::DivideAndConquer()
{
    //  debug timers
    wtime_fragment = microseconds::zero();
    wtime_total = microseconds::zero();
}
DivideAndConquer::~DivideAndConquer()
{
    //  print debug timers
    if (true)
    {
        printf("DaC Fragment Time: %.6f \n", duration_cast<microseconds>(wtime_fragment).count() / 1e6 );
        printf("DaC Total Time:    %.6f \n", duration_cast<microseconds>(wtime_total).count() / 1e6  );
    }
}

void DivideAndConquer::solve(vec_d &Coulomb_mat, vec_d &pot_vec, std::vector<vec_i> &fragments_in, vec_d &delta_rho, vec_d &exponents)
{

    steady_clock::time_point wtime_fragment_start;
    steady_clock::time_point wtime_total_start = steady_clock::now();

    assign_fragments(fragments_in);

    int dim = pot_vec.size();
    vec_d pot_vec_all(dim, 0.0);
    // vec_d delta_rho(dim, 0.0);
    vec_d delta_rho_old(dim, 0.0);
    vec_d lambdas(fragments.size(), 0.0);
    vec_d lambdas_full(dim);
    if (delta_rho.size() != dim)
        delta_rho.resize(dim, 0.0);
    double prev_energy = 1e20;
   

    for(int round_n = 0; round_n < 20; round_n++)
    {
        //  store interaction with all other sites: pot_vec_all = Coulomb_mat @ delta_rho
        cblas_dsymv(CblasRowMajor, CblasUpper, dim, 1.0, &Coulomb_mat[0], dim, &delta_rho[0], 1, 0.0, &pot_vec_all[0], 1);

        wtime_fragment_start = steady_clock::now();
        //#pragma omp parallel for num_threads(1)
        //#pragma omp parallel for
        for(int n = 0; n < (int)fragments.size(); n++)
        {
            int num_sites_n = (int)fragments[n].size();
            std::vector<vec_d> constraints;
            vec_d constraint_vals;
            constraints.push_back(vec_d(num_sites_n, 1.0));
            constraint_vals.push_back(0.0);

            //  solve the minimization problem for this fragment
            //  and keep trying until abs. of all populations are
            //  less than some tollerance
                int num_constr = (int)constraints.size();
                int dim_n = num_sites_n + num_constr;
                vec_d A_mat(dim_n*dim_n, 0.0);
                vec_d b_vec(dim_n);
                vec_d delta_rho_lam_n(dim_n, 0.0);
                vec_d rhs(dim_n);

                if (true)
                {
                    //  copy over Coulomb block and external potential 
                    //  corresponding to this fragment
                    for(int i = 0; i < num_sites_n; i++)
                    {
                        int idx_i = fragments[n][i];
                        for(int j = 0; j < num_sites_n; j++)
                        {
                            int idx_j = fragments[n][j];
                            int orig_idx = idx_i*dim + idx_j;
                            A_mat[i*dim_n + j] = Coulomb_mat[orig_idx];
                        }
                        b_vec[i] = pot_vec[idx_i] + pot_vec_all[idx_i];
                        delta_rho_lam_n[i] = delta_rho[idx_i];
                    }

                    //  copy over constraint as rows
                    for(int i = 0; i < num_constr; i++)
                    {
                        std::copy(constraints[i].begin(), constraints[i].end(), A_mat.begin() + dim_n*num_sites_n + i*dim_n);
                        rhs[num_sites_n + i] = constraint_vals[i];
                    }

                    //  copy over constraint as columns
                    for(int i = 0; i < num_sites_n; i++)
                        for(int j = 0; j < num_constr; j++)
                            A_mat[dim_n*i + num_sites_n + j] = constraints[j][i];
                }

                else
                {
                    //  copy over block and external potential 
                    //  corresponding to this fragment
                    for(int i = 0; i < dim_n; i++)
                    {
                        int idx_i = fragments[n][i];
                        for(int j = 0; j < dim_n; j++)
                        {
                            if ((i == num_sites_n) || (j == num_sites_n))
                                A_mat[i*dim_n + j] = 1.0;
                            else
                            {
                                int idx_j = fragments[n][j];
                                int orig_idx = idx_i*dim + idx_j;
                                A_mat[i*dim_n + j] = Coulomb_mat[orig_idx];
                            }
                        }
                        if (i < num_sites_n)
                        {
                            b_vec[i] = pot_vec[idx_i] + pot_vec_all[idx_i];
                            delta_rho_lam_n[i] = delta_rho[idx_i];
                        }
                    }
                    A_mat[num_sites_n*dim_n + num_sites_n] = 0.0;
                    delta_rho_lam_n[num_sites_n] = 0.0;
                    rhs[num_sites_n] = 0.0;
                }

                //  subtract out the potential from delta_rho on this current fragment
                //  and copy over right hand side (rhs)
                vec_d fragment_pot(dim_n);
                
                //openblas_set_num_threads(1);
                cblas_dsymv(CblasRowMajor, CblasUpper, dim_n, 1.0, &A_mat[0], dim_n, &delta_rho_lam_n[0], 1, 0.0, &fragment_pot[0], 1);
                for(int i = 0; i < num_sites_n; i++)
                {
                    b_vec[i] -= fragment_pot[i];
                    rhs[i] = -b_vec[i];
                }

                //  minimize delta_rho for this fragment
                int ipiv[dim_n];
                int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim_n, 1, &A_mat[0], dim_n, ipiv, &rhs[0], dim_n);

                //  copy over delta_rho and lagrange multipliers to entire collection
                //  also check for maximum population change
                double max_delta_rho = 0.0;
                int max_delta_rho_idx = 0;
                lambdas[n] = rhs[num_sites_n];
                for(int i = 0; i < num_sites_n; i++)
                {
                    delta_rho[fragments[n][i]] = rhs[i];
                    if (abs(rhs[i]) > abs(max_delta_rho))
                    {
                        max_delta_rho = rhs[i];
                        max_delta_rho_idx = i;
                    }
                    lambdas_full[fragments[n][i]] = lambdas[n];
                }

                {
                    vec_d eye(num_sites_n*num_sites_n, 0.0);
                    vec_d x0(num_sites_n);
                    for(int kk = 0; kk < num_sites_n; kk++)
                    {
                        eye[kk*num_sites_n + kk] = 1.0;
                        x0[kk] = rhs[kk];
                    }
                    project_constraints(eye, x0, constraints, constraint_vals, 0.5);
                    for(int i = 0; i < num_sites_n; i++)
                        delta_rho[fragments[n][i]] = x0[i];
                }
        } // end loop over fragments
        wtime_fragment += (steady_clock::now() - wtime_fragment_start);

        //  get rms change in delta_rho
        double sum = 0.0;
        double max_diff = 0.0;
        for(int i = 0; i < dim; i++)
        {
            double diff = delta_rho_old[i] - delta_rho[i];
            sum += diff*diff;
            delta_rho_old[i] = delta_rho[i];
            max_diff = std::max(max_diff, std::abs(diff));
        }
        double rms = sqrt(sum/(double)dim);
        
        //  perform J_mat * rho -> y
        //  perform 0.5 * rho * y = 0.5 * rho * (Jmat * rho)
        //  perform rho * pot
        vec_d y_vec(dim);
        cblas_dsymv(CblasRowMajor, CblasUpper, dim, 1.0, &Coulomb_mat[0], dim, &delta_rho[0], 1, 0.0, &y_vec[0], 1);
        double term1 = 0.5*cblas_ddot(dim, &delta_rho[0], 1, &y_vec[0], 1);
        double term2 = cblas_ddot(dim, &delta_rho[0], 1, &pot_vec[0], 1);
        double energy = term1 + term2;
        double energy_diff = (energy - prev_energy)*2625.5009;
        
        double max_gradient = 0.0;
        for(int i = 0; i < dim; i++)
        {
            double gradient = y_vec[i] + pot_vec[i] + lambdas_full[i];
            max_gradient = std::max(max_gradient, std::abs(gradient));
            //printf("GRADIENT: %d  %.3f  %.3f  %.3f  %.3f\n", i, pot_vec[i]*2625.5009, lambdas_full[i]*2625.5009, y_vec[i]*2625.5009, max_gradient*2625.5009);
        }

        // printf("round %2d: %15.8f %15.8f %15.8f  %15.8f\n", 
        //           round_n, rms, max_diff, energy*2625.5009, max_gradient*2625.5009);
        // std::cin.get();

        prev_energy = energy;
        if (std::abs(energy_diff) < 1.0)
            break;
        if(max_gradient*2625.5009 < 0.5)
            break;
        
    }
    
    //  update total time
    wtime_total += (steady_clock::now() - wtime_total_start);
}

double DivideAndConquer::sign(double x)
{
    return (0 < x) - (x < 0);
}

void DivideAndConquer::project_constraints(vec_d &metric, vec_d &guess_x, std::vector<vec_d> &constraints_in, vec_d &constraint_vals_in, double max_abs_val)
{
    int dim_x = (int)guess_x.size();
    vec_d x0(guess_x);
    double bound = abs(max_abs_val);
    double error = bound*0.0001;
    
    for(int trial = 0; trial < 10; trial++)
    {
        //  constrain any populations over max_abs_val
        std::vector<vec_d> constraints(constraints_in);
        vec_d constraint_vals(constraint_vals_in);
        vec_d orig_x0(x0);
        for(int i = 0; i < dim_x; i++)
        {
            double diff = abs(x0[i]) - bound;
            if (diff > error)
            {
                //printf("FOUND ONE trial=%d:  idx=%d  val=%.3f  diff=%.5e\n", trial, i, x0[i], abs(x0[i]) - abs(max_abs_val));
                vec_d new_constraint(dim_x);
                new_constraint[i] = 1;
                constraints.push_back(new_constraint);
                constraint_vals.push_back(sign(x0[i])*max_abs_val);
                x0.push_back(0.0);
            }
        }
        if(constraint_vals.size() == constraint_vals_in.size())
            break;
            //  no new constraints were added, we are done 

        //  copy over metric
        int dim_full = dim_x + constraint_vals.size();
        vec_d A_mat(dim_full*dim_full);
        for(int i = 0; i < dim_x; i++)
        {
            for(int j = 0; j < dim_x; j++)
                A_mat[i*dim_full + j] = metric[i*dim_x + j];
        }

        //  solve the optimization
        int ipiv2[dim_full];
        apply_constraints(A_mat, x0, constraints, constraint_vals);
        int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim_full, 1, &A_mat[0], dim_full, ipiv2, &x0[0], dim_full);

        // printf("TRIAL:\n");
        // for(int i = 0; i < dim_x; i++)
        //     printf("DELTA: %2d  %8.3f  %8.3f  %8.3f\n", i, orig_x0[i], x0[i], x0[i] - orig_x0[i]);
        // std::cin.get();
    }
    std::copy(&(x0[0]), &(x0[dim_x]), guess_x.begin());
}

void DivideAndConquer::apply_constraints(vec_d &A_mat, vec_d &rhs, std::vector<vec_d> &constraints, vec_d &constraint_vals)
{
    int dim_1 = constraints[0].size();
    int dim_2 = (int)constraints.size();
    int dim_12 = dim_1 + dim_2; // num_sites

    //  copy over constraint as rows
    for(int i = 0; i < dim_2; i++)
    {
        std::copy(constraints[i].begin(), constraints[i].end(), A_mat.begin() + dim_12*dim_1 + i*dim_12);
        rhs[dim_1 + i] = constraint_vals[i];
    }

    //  copy over constraint as columns
    for(int i = 0; i < dim_1; i++)
        for(int j = 0; j < dim_2; j++)
            A_mat[dim_12*i + dim_1 + j] = constraints[j][i];
}

void DivideAndConquer::solve_restrained(vec_d &A_mat_in, vec_d &rhs_in, vec_d &guess, int n_sites, double rho_max, vec_d &out)
{
    /*  Iteratively solve constraint problem. Development only.
    */
    int dim_n = (int)guess.size();
    int ipiv[dim_n];
    double total_cost;
    vec_d rhs(rhs_in);
    vec_d A_mat(A_mat_in);
    vec_d solution(guess);
    for(int n_try = 0; n_try < 10; n_try++)
    {
        std::copy(rhs_in.begin(), rhs_in.end(), rhs.begin());
        std::copy(A_mat_in.begin(), A_mat_in.end(), A_mat.begin());
        total_cost = 0;
        for(int i = 0; i < n_sites; i++)
        {
            double dx = solution[i] - rho_max;
            if(dx > 0)
            {
                printf("OVER MAX: %2d  %8.3f\n", i, solution[i]);
                double cost = 0.04*dx*dx*dx*dx;
                rhs[i] -= cost*4/dx;
                total_cost += cost;
            }
        }

        int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim_n, 1, &A_mat[0], dim_n, ipiv, &rhs[0], dim_n);
        std::copy(rhs.begin(), rhs.end(), solution.begin());

        printf("     ROUND  %2d: %.10f \n", (n_try+1), total_cost);
        for(int i = 0; i < n_sites; i++)
            printf("DELTA2: %2d  %8.3f \n", i, solution[i]);
        std::cin.get();
    }

    std::copy(solution.begin(), solution.end(), out.begin());
}

void DivideAndConquer::assign_fragments(std::vector<vec_i> &fragments_in)
{
    if (fragments.size() != fragments_in.size())
    {
        fragments.clear();
        fragments.resize(fragments_in.size());
        for (size_t i = 0; i < fragments_in.size(); i ++)
        {
            fragments[i].resize(fragments_in[i].size());
            for (size_t j = 0; j < fragments_in[i].size(); j++)
            {
                fragments[i][j] = fragments_in[i][j];
            }
        }
    }
}