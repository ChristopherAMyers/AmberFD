#include "FlucDens.h"
using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::set;

FlucDens::FlucDens()
{
    printf("Inside FlucDens Constructor \n");
}


FlucDens::FlucDens(const int num_sites, 
                   const double *frozen_chg_in, 
                   const double *nuclei_in, 
                   const double *frozen_exp_in, 
                   const double *dynamic_exp_in)
{
    n_sites = num_sites;
    frozen_pop.reserve(num_sites);

    //  copy over parameters
    frozen_chg.assign(frozen_chg_in, frozen_chg_in + n_sites);
    frozen_exp.assign(frozen_exp_in, frozen_exp_in + n_sites);
    dynamic_exp.assign(dynamic_exp_in, dynamic_exp_in + n_sites);
    nuclei.assign(nuclei_in, nuclei_in + n_sites);
    for (size_t i = 0; i < n_sites; i++)
        frozen_pop.push_back(nuclei[i] - frozen_chg[i]);

    //  param_data is for easy access to parameters my names
    param_data["frozen_charges"] = &frozen_chg;
    param_data["frozen_pop"] =     &frozen_pop;
    param_data["frozen_exp"] =     &frozen_exp;
    param_data["dynamic_exp"] =    &dynamic_exp;
    param_data["nuclei"] =         &nuclei;

    //  density distance cutoff info
    dens_cutoff_power_law = dens_cutoff_a*log(dens_cutoff_pct_error) + dens_cutoff_b;

    //  for dynamic densities
    delta_rho_coulomb_mat.resize(n_sites*n_sites, 0.0);
    delta_rho_pot_vec.resize(n_sites, 0.0);
    for (int i = 0; i < n_sites; i++)
        delta_rho_coulomb_mat[i*n_sites + i] = dynamic_exp[i]*5/16;

    //  initialize exclusions
    exclusions_del_frz.resize(num_sites);
    exclusions_frz_frz.resize(num_sites);

    //  initially assign all sites the the same fragment ID
    site_frag_ids.resize(n_sites, -1);
    n_fragments = 0;
    total_time = 0.0;
    use_frag_constraints = true;
    damp_exponent = 1.26697;
    damp_coeff = 0.81479;
    damp_sum.resize(n_sites, 0.0);
}

FlucDens::~FlucDens()
{}

double FlucDens::dot3Vec(const vec_d &coords, int i, int j)
{
    double x = coords[i]     - coords[j];
    double y = coords[i + 1] - coords[j + 1];
    double z = coords[i + 2] - coords[j + 2];
    return x*x + y*y + z*z;
}

void FlucDens::set_frag_constraints(const bool constr_frags)
{
    /* choose to constrain total charge per fragment (True) or as a whole system (False) */
    use_frag_constraints = constr_frags;
}

void FlucDens::create_del_exclusions_from_frgment(const std::vector<int> frag_idx)
{
    for (auto del_i: frag_idx)
        for(auto frz_j: frag_idx)
            add_del_frz_exclusion(del_i, frz_j);
}

void FlucDens::add_del_frz_exclusion(int delta_i, int frz_j)
{
    /*  delta_rho_i and frozen_rho_j interactions
        will be excluded */
    if ((delta_i > n_sites) || (frz_j > n_sites) || (delta_i < 0) || (frz_j < 0))
        throw out_of_bounds_eror("add_del_frz_exclusion()", delta_i, frz_j);
    exclusions_del_frz[delta_i].push_back(frz_j);
}

void FlucDens::get_del_frz_exclusions(const int particle1, std::vector<int> &particles2) const
{
    if (particle1 < 0 || particle1 >= (int) exclusions_frz_frz.size()) 
        throw "Index out of range";
    particles2 = exclusions_del_frz[particle1];
}

void FlucDens::add_frz_frz_exclusion(int frz_i, int frz_j)
{
    /*  frozen_rho_i and frozen_rho_j interactions
        will be excluded */
    if ((frz_i > n_sites) || (frz_j > n_sites) || (frz_i < 0) || (frz_j < 0))
    {   
        char buffer[100];
        sprintf(buffer, "Frozen exception index pair (%d,%d) out of bounds", frz_i, frz_j);
        throw std::out_of_range(buffer);
    }
    exclusions_frz_frz[frz_i].push_back(frz_j);
    exclusions_frz_frz[frz_j].push_back(frz_i);
}

int FlucDens::get_num_frz_frz_exclusions() const
{
    return exclusions_frz_frz.size(); 
}

void FlucDens::get_frz_frz_exclusions(const int particle1, std::vector<int> &particles2) const
{
    if (particle1 < 0 || particle1 >= (int) exclusions_frz_frz.size()) 
        throw "Index out of range";
    particles2 = exclusions_frz_frz[particle1];
}

void FlucDens::create_frz_exclusions_from_bonds(const vector<pair<int, int> > bonds, int bond_cutoff)
{
    /*  determine exclusions for frozen-frozen interactions from bonds up to 
        an arbitrary bond neighbor length. Edited from OpenMM */ 
    if (bond_cutoff < 1)
        return;
    for (auto& bond : bonds)
    {
        if (bond.first < 0 || bond.second < 0 || bond.first >= n_sites || bond.second >= n_sites)
            throw "createExclusionsFromBonds: Illegal particle index in list of bonds";
    }
    vector<std::set<int> > exclusions(n_sites);
    vector<std::set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        int p1 = bond.first;
        int p2 = bond.second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bond_cutoff-1; level++) {
        vector<std::set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) n_sites; i++)
            for (int j : currentExclusions[i])
                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i)
                add_frz_frz_exclusion(i, j);    
}

void FlucDens::add_fragment(const std::vector<int> site_idx_list)
{
    for(auto delta_i: site_idx_list)
    {
        if ((delta_i > n_sites) || (delta_i < 0) )
            throw out_of_bounds_eror("add_fragment()", delta_i);
        
        site_frag_ids[delta_i] = n_fragments;
        for(auto frz_j: site_idx_list)
            add_del_frz_exclusion(delta_i, frz_j);
    }
    n_fragments += 1;
}

void FlucDens::change_dyn_exp(const int index, const double value)
{
    dynamic_exp[index] = value;
    delta_rho_coulomb_mat[index*n_sites + index] = dynamic_exp[index]*5/16;
}

void FlucDens::change_frz_exp(const int index, const double value)
{
    frozen_exp[index] = value;
}

void FlucDens::set_dampening(double coeff, double exponent)
{
    damp_coeff = coeff;
    damp_exponent = exponent;
}

double FlucDens::get_total_time()
{
    return total_time;
}

double FlucDens::calc_overlap(const vec_d &coords)
{
    size_t i, j;
    double overlap = 0.0;
    double r, r2, inv_r;

    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            if ((nuclei[i] < 1) || (nuclei[j] < 1))
                continue;
                
            //  symmetric Coulomb matrix elements in row_major format
            size_t matrix_idx_1 = i*n_sites + j;
            size_t matrix_idx_2 = j*n_sites + i;

            //  distances and inverse distances
            r2 = dot3Vec(coords, (int)j*3, (int)i*3);
            r = sqrt(r2);
            inv_r = 1/r;

            //  extract density exponentials
            double a = frozen_exp[i];
            double b = frozen_exp[j];

            if (std::find(exclusions_frz_frz[i].begin(), exclusions_frz_frz[i].end(), j) == exclusions_frz_frz[i].end())
            {
                if (use_long_range_approx(r, a, b))
                {            }
                else if (true)
                {
                    //  compute exponentials only once
                    double exp_ar = exp(-a*r);
                    double exp_br = exp(-b*r);

                    overlap += frz_frz_overlap(inv_r, a, b, exp_ar, exp_br)*frozen_pop[i]*frozen_pop[j];
                }
            }
        }
    }

    return overlap;
}

double FlucDens::frz_frz_overlap(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br)
{
    /* Calculate the overlap between two atoms */

    double ab_diff = b - a;
    if (abs(ab_diff) > 0.001)
    {
        double a2 = a*a;
        double b2 = b*b;
        double ab4 = a2*a2*b2*b2;
        double r = 1/inv_r;
        double denom = a2 - b2;
        double denom2 = denom*denom;

        double c1 = 4*ab4/(denom*denom2)*inv_r;
        double c2 = ab4/a/denom2;
        double c3 = -c1;
        double c4 = c2*a/b;

        return (exp_ar*(c1 + c2) + exp_br*(c3 + c4))/(8*M_PI);
    }
    else
    {
        double ab_diff = b - a;
        double ar = a/inv_r;
        double ar2 = ar*ar;
        double ar3 = ar2*ar;

        double term0 = (a*a*(3 + 3*ar + ar2))/192;
        double term1 = (a*ab_diff*(9 + 9*ar + 2*ar2 - ar3))/384;
        double term2 = (ab_diff*ab_diff*ar2*(-5 - 5*ar + ar2))/1280;

        return (term0 + term1 + term2)*a*exp_ar/M_PI;
    }
}

double FlucDens::calc_energy(const vec_d &positions, bool calc_frz, bool calc_pol)
{
    size_t i, j;

    frozen_energy = 0.0;
    polarization_energy = 0.0;

    if (!calc_pol && !calc_frz)
        return 0.0;
    
    std::fill(delta_rho_pot_vec.begin(), delta_rho_pot_vec.end(), 0.0);
    std::fill(damp_sum.begin(), damp_sum.end(), 0.0);

    //#pragma omp parallel for private (frozen_energy)
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances and inverse distances
            double deltaR[Nonbonded::RMaxIdx];
            Nonbonded::calc_dR(positions, (int)j*3, (int)i*3, deltaR);

            calc_one_electro(deltaR, i, j, calc_pol, calc_frz);
        }
    }
    if (calc_pol)
        solve_minimization();

    return frozen_energy + polarization_energy;
}

void FlucDens::calc_one_electro(double* deltaR, int i, int j, bool calc_pol, bool calc_frz)
{
    const double b_frz = frozen_exp[j];
    const double b_del = dynamic_exp[j];
    const double a_frz = frozen_exp[i];
    const double a_del = dynamic_exp[i];
    
    //  extract deltaR
    const double r = deltaR[Nonbonded::RIdx];
    const double inv_r = deltaR[Nonbonded::RInvIdx];

    //  symmetric Coulomb matrix elements
    const size_t matrix_idx_1 = i*n_sites + j;
    const size_t matrix_idx_2 = j*n_sites + i;

    if (use_long_range_approx(r, a_frz, b_frz)
    && use_long_range_approx(r, a_del, b_del) )
    {            }
    else if (true)
    {
        const double exp_ar_del = exp(-a_del*r);
        const double exp_ar_frz = exp(-a_frz*r);
        const double exp_br_frz = exp(-b_frz*r);
        const double exp_br_del = exp(-b_del*r);

        if (calc_pol)
        {
            //  delta_rho - delta_rho
            const double del_del_ee = elec_elec_penetration(inv_r, a_del, b_del, exp_ar_del, exp_br_del);
            delta_rho_coulomb_mat[matrix_idx_1] = del_del_ee;
            delta_rho_coulomb_mat[matrix_idx_2] = del_del_ee;

            //  delta_rho - frozen_rho
            if (std::find(exclusions_del_frz[i].begin(), exclusions_del_frz[i].end(), j) == exclusions_del_frz[i].end())
            {
                
                const double del_frz_ee  = elec_elec_penetration(inv_r, a_del, b_frz, exp_ar_del, exp_br_frz);
                const double frz_del_ee  = elec_elec_penetration(inv_r, a_frz, b_del, exp_ar_frz, exp_br_del);
                const double del_a_nuc_b = elec_nuclei_pen(inv_r, a_del, exp_ar_del);
                const double del_b_nuc_a = elec_nuclei_pen(inv_r, b_del, exp_br_del);
                
                delta_rho_pot_vec[i] += frozen_pop[j]*del_frz_ee - nuclei[j]*del_a_nuc_b;
                delta_rho_pot_vec[j] += frozen_pop[i]*frz_del_ee - nuclei[i]*del_b_nuc_a;

                damp_sum[i] += exp(-damp_exponent*r);
                damp_sum[j] += exp(-damp_exponent*r);
            }
        }
        //  frozen_rho - frozen_rho
        if (calc_frz)
        if (std::find(exclusions_frz_frz[i].begin(), exclusions_frz_frz[i].end(), j) == exclusions_frz_frz[i].end())
        {
            const double frz_frz_ee  = elec_elec_penetration(inv_r, a_frz, b_frz, exp_ar_frz, exp_br_frz);
            const double frz_a_nuc_b = elec_nuclei_pen(inv_r, a_frz, exp_ar_frz);
            const double frz_b_nuc_a = elec_nuclei_pen(inv_r, b_frz, exp_br_frz);
            const double nuc_nuc = nuclei[i]*nuclei[j]*inv_r;

            frozen_energy +=  frozen_pop[i]*frozen_pop[j]*frz_frz_ee 
                            - frozen_pop[i]*nuclei[j]*frz_a_nuc_b
                            - frozen_pop[j]*nuclei[i]*frz_b_nuc_a
                            + nuc_nuc;
        }
    }
}

double FlucDens::elec_nuclei_pen(const double inv_r, const double a, const double exp_ar)
{
    double nuc_elec;
    if (inv_r > 1E4)
        nuc_elec = a*0.5;
    else
    {
        nuc_elec = inv_r - exp_ar*(inv_r + a*0.5);
    }
    return nuc_elec;
}


double FlucDens::elec_elec_penetration(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br)
{
    /*  electron electron charge penetration term */
    double ab_diff = b - a;
    if (abs(ab_diff) > 0.001)
    {
        double a2 = a*a; double a4 = a2*a2; double a6 = a2*a4;
        double b2 = b*b; double b4 = b2*b2; double b6 = b2*b4;
        double denom1 = a2 - b2;
        double denom2 = -denom1;
        double den1_2 = denom1*denom1;
        double den2_2 = denom2*denom2;
        double c1 = a*b4 / (2*den1_2);
        double c2 = (b6 - 3*a2*b4) / (den1_2*denom1);
        double c3 = b*a4 / (2*den2_2);
        double c4 = (a6 - 3*b2*a4) / (den2_2*denom2);

        return inv_r - exp_ar*(c1 - c2*inv_r) - exp_br*(c3 - c4*inv_r);
    }
    else
    {
        double ar = a/inv_r;
        double ar2 = ar*ar;
        double ar3 = ar2*ar;
        double r = 1/inv_r;
        double term0 = (-48 - 33*ar - 9*ar2 - ar3)/(48.*r);
        double term1 = (15 + 15*ar + 6*ar2 + ar3)/96;
        double term2 = -((30 + 30*ar + 15*ar2 + 5*ar3 + ar2*ar2)*r)/(320.*ar);

        return inv_r + exp_ar*(term0 + term1*ab_diff + term2*ab_diff*ab_diff);
    }
}


void FlucDens::assign_constraints()
{
    if (use_frag_constraints)
    {
        if ((int)constraints.size() != n_fragments)
        {
            constraints.resize(n_fragments, std::vector<int>(n_sites, 0));
            for(int i = 0; i < n_sites; i ++)
            {
                int frag_id = site_frag_ids[i];
                if (frag_id == -1)
                    throw std::runtime_error("assign_constraints: Site " + std::to_string(i+1) + " not assigned a fragment.");
                constraints[frag_id][i] = 1;
            }
        }
    }
    else
    {
        if ((int)constraints.size() != 1)
        {
            //printf("RESIZING CONSTRAINTS %d  %d\n", (int)constraints.size(), (int)use_frag_constraints);
            constraints.clear();
            constraints.resize(1, std::vector<int>(n_sites, 1));
        }
    }
}

void FlucDens::solve_minimization()
{
    assign_constraints();
    int n_constr = (int)constraints.size();
    int dim = n_sites + n_constr;
    int info;
    vec_d A_mat(dim*dim, 0.0);
    vec_d B_vec(dim, 0);

    //  copy over dampening
    for (int i = 0; i < n_sites; i++)
        delta_rho_coulomb_mat[i*n_sites + i] = dynamic_exp[i]*5/16 + damp_sum[i]*damp_coeff;
    
    //  copy over potential vector
    vec_d neg_pot = delta_rho_pot_vec;
    cblas_dscal(n_sites, -1, &neg_pot[0], 1);   //negate potential (right hand side of matrix equation)
    std::copy(neg_pot.begin(), neg_pot.end(), B_vec.begin());

    //  copy over Coulomb matrix
    vec_d::iterator start = delta_rho_coulomb_mat.begin();
    for(int i = 0; i < n_sites; i++)
        std::copy(start + i*n_sites, start + (i+1)*n_sites, A_mat.begin() + i*dim);
    
    //  copy over constraint as rows
    for(int i = 0; i < n_constr; i++)
        std::copy(constraints[i].begin(), constraints[i].end(), A_mat.begin() + dim*n_sites + i*dim);

    //  copy over constraint as columns
    for(int i = 0; i < n_sites; i++)
        for(int j = 0; j < n_constr; j++)
            A_mat[dim*i + n_sites + j] = constraints[j][i];

    A_mat_save = A_mat;
    B_vec_save = B_vec;

    //  solve matrix equation A_mat * x = B_vec for x
    int ipiv[dim];
    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim, 1, &A_mat[0], dim, ipiv, &B_vec[0], dim);
    delta_rho.assign(B_vec.begin(), B_vec.begin() + n_sites);

    //  Since all constraints are linear in delta_rho and result in sums equal to zero,
    //  the polarization energy is simply 0.5 * delta_rho^T * pot_vec
    polarization_energy = 0.5*cblas_ddot(n_sites, &delta_rho[0], 1, &delta_rho_pot_vec[0], 1);
    return;

    /*
    //  solve for polarization energy = 0.5 * dP^T Coul_mat * dP + pot_vec * dP
    //  energy is simplified down to  = -0.5 * dP^T (constr_mat^T * lambda + neg_pot)
    vec_d constraint_flat(constraints.size()*n_sites, 0.0); // in column-major format
    vec_d lamb(B_vec.begin() + n_sites, B_vec.end());
    for(int i = 0; i < n_constr; i++)
        std::copy(constraints[i].begin(), constraints[i].end(), constraint_flat.begin() + i*n_sites);
    
    //  The constr_mat^T * lambda + neg_pot part; result stored in neg_pot
    cblas_dgemv(CblasColMajor, CblasNoTrans, n_sites, n_constr, 1, &constraint_flat[0], n_sites, &lamb[0], 1, 1, &neg_pot[0], 1);
    //  -0.5 * dP^T (previous_result) part
    polarization_energy = -0.5*cblas_ddot(n_sites, &delta_rho[0], 1, &neg_pot[0], 1);
    */
}

bool FlucDens::use_long_range_approx(double r, double a, double b)
{
    return false;
}

vec_d FlucDens::get_params_by_name(const std::string param_name)
{
    return vec_d(param_data[param_name]->begin(), param_data[param_name]->end());
}

void FlucDens::print_params(const std::string message, const std::string param_name)
{
    vec_d *params = param_data[param_name];
    size_t i = 0;
    int n_per_row = std::min(n_sites, (size_t)5);
    printf("\n %s\n", message.c_str());
    printf(" Number of elements: %d \n", (int)frozen_exp.size());
    for (i = 0; i < n_sites; i ++)
    {
        printf(" %10.5f", (*params)[i]);
        if (i % n_per_row == (n_per_row - 1))
            printf("\n");
    }
    printf("\n");
}

vec_d FlucDens::get_rho_coulomb_mat()
{
    return vec_d(delta_rho_coulomb_mat);
}

vec_d FlucDens::get_rho_pot_vec()
{
    return vec_d(delta_rho_pot_vec);
}

vec_d FlucDens::get_delta_rho()
{
    return delta_rho;
}

double FlucDens::get_frozen_energy()
{
    return frozen_energy;
}

double FlucDens::get_polarization_energy()
{
    return polarization_energy;
}

std::vector<std::vector<int>> FlucDens::get_constraints()
{
    return constraints;
}

std::out_of_range FlucDens::out_of_bounds_eror(const char *msg, const int idx1)
{
    char buffer[100];
    sprintf(buffer, "Index %d out of bounds: %s", idx1, msg);
    return std::out_of_range(buffer);
}

std::out_of_range FlucDens::out_of_bounds_eror(const char *msg, const int idx1, const int idx2)
{
    char buffer[100];
    sprintf(buffer, "Index pair (%d,%d) out of bounds: %s", idx1, idx2, msg);
    return std::out_of_range(buffer);
}