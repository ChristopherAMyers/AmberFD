#include "LinearSolvers.h"

DivideAndConquer::DivideAndConquer()
{
    
}

void DivideAndConquer::solve(vec_d &Coulomb_mat, vec_d &pot_vec, std::vector<vec_i> &fragments_in)
{
    assign_fragments(fragments_in);

    int dim = pot_vec.size();
    vec_d pot_vec_all(dim, 0.0);
    vec_d delta_rho(dim, 0.0);
    vec_d delta_rho_old(dim, 0.0);

    for(int round_n = 0; round_n < 15; round_n++)
    {
        //printf("Round 1\n");
        //  store interaction with all other sites: pot_vec_all = Coulomb_mat @ delta_rho
        cblas_dsymv(CblasRowMajor, CblasUpper, dim, 1.0, &Coulomb_mat[0], dim, &delta_rho[0], 1, 0.0, &pot_vec_all[0], 1);
        for(int n = 0; n < (int)fragments.size(); n++)
        {
            //printf("Fragment %d\n", n);
            int dim_n = (int)fragments[n].size();

            vec_d A_mat(dim_n*dim_n);
            vec_d b_vec(dim_n);
            vec_d delta_rho_n(dim_n);

            //  copy over block and external potential 
            //  corresponding to this fragment
            for(int i = 0; i < dim_n; i++)
            {
                int idx_i = fragments[n][i];
                for(int j = 0; j < dim_n; j++)
                {
                    int idx_j = fragments[n][j];
                    int orig_idx = idx_i*dim + idx_j;
                    A_mat[i*dim_n + j] = Coulomb_mat[orig_idx];
                }
                b_vec[i] = pot_vec[idx_i] + pot_vec_all[idx_i];
                delta_rho_n[i] = delta_rho[idx_i];
            }

            //  subtract out the potential from delta_rho on this current fragment
            //  and copy over right hand side (rhs)
            vec_d fragment_pot(dim_n);
            vec_d rhs(dim_n);
            cblas_dsymv(CblasRowMajor, CblasUpper, dim_n, 1.0, &A_mat[0], dim_n, &delta_rho_n[0], 1, 0.0, &fragment_pot[0], 1);
            for(int i = 0; i < dim_n; i++)
            {
                b_vec[i] -= fragment_pot[i];
                rhs[i] = -b_vec[i];
            }

            //  minimize delta_rho for this fragment
            int ipiv[dim_n];
            int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim_n, 1, &A_mat[0], dim_n, ipiv, &rhs[0], dim_n);

            //  copy over delta_rho to entire collection
            for(int i = 0; i < dim_n; i++)
                delta_rho[fragments[n][i]] = rhs[i];
        }

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

        printf("round %2d: %15.8f %15.8f %15.8f\n", round_n, rms, max_diff, energy*2625.5009);
    }
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