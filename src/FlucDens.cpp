#include "FlucDens.h"
using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::set;

FlucDens::FlucDens(const int num_sites, 
                   const double *frozen_chg_in, 
                   const double *nuclei_in, 
                   const double *frozen_exp_in, 
                   const double *dynamic_exp_in)
{
    n_sites = num_sites;
    frozen_pop.reserve(num_sites);
    frozen_forces.resize(num_sites);
    total_forces.resize(num_sites);

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
    J_mat.resize(n_sites*n_sites, 0.0);
    // dJ_dx.resize(n_sites*n_sites, 0.0);
    // dJ_dy.resize(n_sites*n_sites, 0.0);
    // dJ_dz.resize(n_sites*n_sites, 0.0);
    dJ_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dPot_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dPot_dPos_trans.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dDamp_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));

    pot_vec.resize(n_sites, 0.0);
    for (int i = 0; i < n_sites; i++)
        J_mat[i*n_sites + i] = dynamic_exp[i]*5/16;

    //  initialize exclusions
    exclusions_del_frz.resize(num_sites);
    exclusions_frz_frz.resize(num_sites);

    //  initially assign all sites the the same fragment ID
    site_frag_ids.resize(n_sites, -1);
    n_fragments = 0;
    total_time = 0.0;
    use_frag_constraints = true;
    damp_exponent = 1.0;
    damp_coeff = 0.0;
    damp_sum.resize(n_sites, 0.0);
    calc_forces = true;

    //  initialize energies
    initialize_calculation();
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

void FlucDens::create_del_exclusions_from_fragment(const std::vector<int> frag_idx)
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
    exclusions_del_frz[delta_i].insert(frz_j);
}

std::set<int> FlucDens::get_del_frz_exclusions(const int particle1) const
{
    std::set<int> particles2;
    if (particle1 < 0 || particle1 >= (int) exclusions_frz_frz.size()) 
        throw "Index out of range";
    particles2 = exclusions_del_frz[particle1];
    return particles2;
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
    exclusions_frz_frz[frz_i].insert(frz_j);
    exclusions_frz_frz[frz_j].insert(frz_i);
}

int FlucDens::get_num_frz_frz_exclusions() const
{
    return exclusions_frz_frz.size(); 
}

std::vector<int> FlucDens::get_frz_frz_exclusions(const int particle1) const
{
    if (particle1 < 0 || particle1 >= (int) exclusions_frz_frz.size()) 
        throw "Index out of range";
    std::vector<int> particles2;
    particles2.resize(exclusions_frz_frz[particle1].size());
    std::copy(exclusions_frz_frz[particle1].begin(), exclusions_frz_frz[particle1].end(), particles2.begin());

    int i = 0;
    for(auto &x: exclusions_frz_frz[particle1])
    {
        printf("%d  %d  %d \n", (int)i, particles2[i], x);
        i++;
    }

    return particles2;
    //particles2 = exclusions_frz_frz[particle1];
}

void FlucDens::create_frz_exclusions_from_bonds(const vector<pair<int, int> > bonds, int bond_cutoff)
{
    vector<std::set<int> > exclusions = Nonbonded::calc_exclusions_from_bonds(bonds, bond_cutoff, n_sites);
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

void FlucDens::set_dyn_exp(const int index, const double value)
{
    dynamic_exp[index] = value;
    J_mat[index*n_sites + index] = dynamic_exp[index]*5/16;
}

void FlucDens::set_dyn_exp(vec_d exponents)
{
    if ((int)exponents.size() != n_sites)
        throw std::runtime_error("radi_list length does not equal n_sites");
    dynamic_exp.assign(exponents.begin(), exponents.end());
}

void FlucDens::set_frz_exp(const int index, const double value)
{
    frozen_exp[index] = value;
}

void FlucDens::set_dampening(double coeff, double exponent)
{
    damp_coeff = coeff;
    damp_exponent = exponent;
}

void FlucDens::set_calc_forces(bool calculate_forces)
{
    calc_forces = calculate_forces;
}

double FlucDens::get_total_time()
{
    return total_time;
}

int FlucDens::get_num_constraints()
{
    return constraints.size();
}

std::vector<vec_d> FlucDens::get_forces()
{
    std::vector<vec_d> rtn(n_sites, vec_d(3, 0.0));
    for(size_t i = 0; i < total_forces.size(); i++)
    {
        rtn[i][0] = total_forces[i][0];
        rtn[i][1] = total_forces[i][1];
        rtn[i][2] = total_forces[i][2];
    }
    return rtn;
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

void FlucDens::initialize_calculation()
{ 
    using std::fill;
    fill(damp_sum.begin(), damp_sum.end(), 0.0);
    fill(pot_vec.begin(), pot_vec.end(), 0.0);
    fill(dJ_dPos.begin(), dJ_dPos.end(), Vec3(0, 0, 0));
    fill(dPot_dPos.begin(), dPot_dPos.end(), Vec3(0, 0, 0));
    fill(dPot_dPos_trans.begin(), dPot_dPos_trans.end(), Vec3(0, 0, 0));
    fill(dDamp_dPos.begin(), dDamp_dPos.end(), Vec3(0, 0, 0));
    fill(frozen_forces.begin(), frozen_forces.end(), Vec3(0, 0, 0));
    fill(total_forces.begin(), total_forces.end(), Vec3(0, 0, 0));
    total_energies.reset();
}

double FlucDens::calc_energy(const vec_d &positions, bool calc_frz, bool calc_pol)
{
    if (!calc_pol && !calc_frz)
        return 0.0;

    initialize_calculation();
    size_t i, j;
    Energies energies;
    energies.reset_all();

    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances and inverse distances
            DeltaR deltaR(positions, (int)i*3, (int)j*3);

            calc_one_electro(deltaR, i, j, calc_pol, calc_frz, energies);
            total_energies.frz += energies.frz;
            total_energies.nuc_nuc += energies.nuc_nuc;
            total_energies.elec_elec += energies.elec_elec;
            total_energies.elec_nuc += energies.elec_nuc;
        }
    }
    if (calc_pol)
        solve_minimization();

    return total_energies.total();
}

void FlucDens::calc_one_electro(DeltaR &deltaR, int i, int j, bool calc_pol, bool calc_frz, Energies& energies)
{
    /*  calculate a single electrostatic interaction between a pair of atoms */
    
    double frozen_energy = 0.0;
    double E_ee = 0.0, E_eZ1 = 0.0, E_eZ2 = 0.0, E_ZZ = 0.0;
    

    const double b_frz = frozen_exp[j];
    const double b_del = dynamic_exp[j];
    const double a_frz = frozen_exp[i];
    const double a_del = dynamic_exp[i];

    //  extract deltaR
    const double r = deltaR.r;
    const double inv_r = deltaR.r_inv;

    //  symmetric Coulomb matrix elements
    const int idx_ij = i*n_sites + j;
    const int idx_ji = j*n_sites + i;

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
            Vec3 dR = deltaR.dR*inv_r;
            double ee_ddR;
            const double del_del_ee = elec_elec_energy(inv_r, a_del, b_del, exp_ar_del, exp_br_del, ee_ddR);
            J_mat[idx_ij] = del_del_ee;
            J_mat[idx_ji] = del_del_ee;
            dJ_dPos[idx_ij] =  ee_ddR * dR;
            dJ_dPos[idx_ji] = -ee_ddR * dR;

            //  delta_rho - frozen_rho
            if (std::find(exclusions_del_frz[i].begin(), exclusions_del_frz[i].end(), j) == exclusions_del_frz[i].end())
            {
                double del_frz_ee_ddR, frz_del_ee_ddR, del_a_nuc_b_ddR, del_b_nuc_a_ddR;
                const double del_frz_ee  = elec_elec_energy(inv_r, a_del, b_frz, exp_ar_del, exp_br_frz, del_frz_ee_ddR);
                const double frz_del_ee  = elec_elec_energy(inv_r, a_frz, b_del, exp_ar_frz, exp_br_del, frz_del_ee_ddR);
                const double del_a_nuc_b = elec_nuclei_energy(inv_r, a_del, exp_ar_del, del_a_nuc_b_ddR);
                const double del_b_nuc_a = elec_nuclei_energy(inv_r, b_del, exp_br_del, del_b_nuc_a_ddR);
                
                pot_vec[i] += frozen_pop[j]*del_frz_ee + nuclei[j]*del_a_nuc_b;
                pot_vec[j] += frozen_pop[i]*frz_del_ee + nuclei[i]*del_b_nuc_a;

                double dampening = exp(-damp_exponent*r);
                damp_sum[i] += dampening;
                damp_sum[j] += dampening;

                double dDamp_dR = -damp_exponent*damp_coeff*dampening;
                dDamp_dPos[idx_ij] =  dDamp_dR*dR;
                dDamp_dPos[idx_ji] = -dDamp_dR*dR;

                
                double dP_dR_ij = frozen_pop[j]*del_frz_ee_ddR + nuclei[j]*del_a_nuc_b_ddR;
                double dP_dR_ji = frozen_pop[i]*frz_del_ee_ddR + nuclei[i]*del_b_nuc_a_ddR;
                dPot_dPos[idx_ij] = dP_dR_ij*dR;
                dPot_dPos[idx_ji] = -dP_dR_ji*dR;
                dPot_dPos_trans[idx_ji] = dP_dR_ij*dR;
                dPot_dPos_trans[idx_ij] = -dP_dR_ji*dR;
            }
        }
        //  frozen_rho - frozen_rho
        if (calc_frz)
        if (std::find(exclusions_frz_frz[i].begin(), exclusions_frz_frz[i].end(), j) == exclusions_frz_frz[i].end())
        {
            double frz_frz_ee_ddR, frz_a_nuc_b_ddR, frz_b_nuc_a_ddR;
            const double frz_frz_ee  = elec_elec_energy(inv_r, a_frz, b_frz, exp_ar_frz, exp_br_frz, frz_frz_ee_ddR);
            const double frz_a_nuc_b = elec_nuclei_energy(inv_r, a_frz, exp_ar_frz, frz_a_nuc_b_ddR);
            const double frz_b_nuc_a = elec_nuclei_energy(inv_r, b_frz, exp_br_frz, frz_b_nuc_a_ddR);
            E_ZZ = nuclei[i]*nuclei[j]*inv_r;
            E_ee = frozen_pop[i]*frozen_pop[j]*frz_frz_ee;
            E_eZ1 = frozen_pop[i]*nuclei[j]*frz_a_nuc_b;
            E_eZ2 = frozen_pop[j]*nuclei[i]*frz_b_nuc_a;

            Vec3 dR_norm = deltaR.dR*inv_r;
            Vec3 d_ee_dR = frozen_pop[i]*frozen_pop[j]*frz_frz_ee_ddR*dR_norm;
            Vec3 d_eZ_dR1 = frozen_pop[i]*nuclei[j]*frz_a_nuc_b_ddR*dR_norm;
            Vec3 d_eZ_dR2 = frozen_pop[j]*nuclei[i]*frz_b_nuc_a_ddR*dR_norm;
            Vec3 d_ZZ_dR = -nuclei[i]*nuclei[j]*inv_r*inv_r*dR_norm;

            Vec3 force_i = -(d_ee_dR + d_eZ_dR1 + d_eZ_dR2 + d_ZZ_dR);
            frozen_forces[i] += force_i;
            frozen_forces[j] -= force_i;

            frozen_energy =  (E_ee + E_eZ1) + (E_ZZ + E_eZ2);
        }
    }
    energies.elec_elec = E_ee;
    energies.elec_nuc = E_eZ1 + E_eZ2;
    energies.nuc_nuc = E_ZZ;
    energies.frz = frozen_energy;
}

double FlucDens::elec_nuclei_energy(const double inv_r, const double a, const double exp_ar, double &dEdR)
{
    double nuc_elec;
    if (inv_r > 1E4)
    {
        nuc_elec = -a*0.5;
        dEdR = 0.0;
    }
    else
    {
        nuc_elec = -inv_r + exp_ar*(inv_r + a*0.5);
        dEdR = -nuc_elec*inv_r - 0.5*a*exp_ar*(inv_r + a);
    }
    return nuc_elec;
}


double FlucDens::elec_elec_energy(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br, double &dEdR)
{
    double E_ee = 0.0;
    /*  electron electron charge penetration term */
    double d = b - a;
    if (abs(d) > 0.001)
    {
        double a2 = a*a; double a4 = a2*a2; double a6 = a2*a4;
        double b2 = b*b; double b4 = b2*b2; double b6 = b2*b4;
        double denom1 = a2 - b2;
        double denom2 = -denom1;
        double den1_2 = denom1*denom1;
        double den2_2 = denom2*denom2;
        double c1 = a*b4 / (2*den1_2);
        double c2 = -(b6 - 3*a2*b4) / (den1_2*denom1);
        double c3 = b*a4 / (2*den2_2);
        double c4 = -(a6 - 3*b2*a4) / (den2_2*denom2);

        E_ee = inv_r - exp_ar*(c1 + c2*inv_r) - exp_br*(c3 + c4*inv_r);
        dEdR = -E_ee*inv_r + exp_ar*(a*c1 + inv_r*(a*c2 - c1)) + exp_br*(b*c3 + inv_r*(b*c4 - c3));
    }
    else
    {
        double r = 1/inv_r;
        double ad = a*d;
        double a2 = a*a;
        double a2_minus_3d = 2*a - 3*d;
        double A0 = (5*ad - 3*d*d - 22*a2)/(32*a);
        double A1 = (5*ad - 3*d*d - 6*a2)/32;
        double A2 = -a*a2_minus_3d*a2_minus_3d/192;
        double A3 = a*ad*a2_minus_3d/192;
        double A4 = -a2*ad*d/320;
        double r2 = r*r;
        double r3 = r*r2;
        double r4 = r2*r2;
        double one_minus_exp = 1 - exp_ar;

        E_ee = inv_r*one_minus_exp + exp_ar*(A0 + r*A1 + r2*A2 + r3*A3 + r4*A4);
        dEdR = -a*E_ee + a*inv_r - inv_r*inv_r*one_minus_exp + exp_ar*(A1 + 2*A2*r + 3*A3*r2 + 4*A4*r3);
    }
    return E_ee;
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
        J_mat[i*n_sites + i] = dynamic_exp[i]*5/16 + damp_sum[i]*damp_coeff;
    
    //  copy over potential vector
    vec_d neg_pot = pot_vec;
    cblas_dscal(n_sites, -1, &neg_pot[0], 1);   //negate potential (right hand side of matrix equation)
    std::copy(neg_pot.begin(), neg_pot.end(), B_vec.begin());

    //  copy over Coulomb matrix
    vec_d::iterator start = J_mat.begin();
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

    //  simplified version of polarization energy
    if (false)
    {
        //  Since all constraints are linear in delta_rho and result in sums equal to zero,
        //  the polarization energy is simply 0.5 * delta_rho^T * pot_vec
        //total_pol_energy = 0.5*cblas_ddot(n_sites, &delta_rho[0], 1, &delta_rho_pot_vec[0], 1);
        total_energies.pol = 0.5*cblas_ddot(n_sites, &delta_rho[0], 1, &pot_vec[0], 1);
    }
    //  extended version, used for debugging parts of polarization energy and forces
    else
    {
        //printf("\n ####   Using extended POL energy   #### \n");
        //  create a temp vector to store result in
        vec_d y_vec(n_sites);

        //  perform J_mat * rho -> y
        cblas_dsymv(CblasRowMajor, CblasUpper, n_sites, 1.0, &J_mat[0], n_sites, &delta_rho[0], 1, 0.0, &y_vec[0], 1);

        //  perform 0.5 * rho * y = 0.5 * rho * (Jmat * rho)
        double term1 = 0.5*cblas_ddot(n_sites, &delta_rho[0], 1, &y_vec[0], 1);

        //  perform rho * pot
        double term2 = cblas_ddot(n_sites, &delta_rho[0], 1, &pot_vec[0], 1);

        total_energies.pol = term1 + term2;
    }

    for(int k = 0; k < n_sites; k++)
    {
        total_forces[k] = Vec3(0, 0, 0);
        for(int j = 0; j < n_sites; j++)
        {
            int kj_idx = k*n_sites + j;
            total_forces[k] -= delta_rho[k]*dJ_dPos[kj_idx]*delta_rho[j];
            total_forces[k] -= (delta_rho[k]*dPot_dPos[kj_idx] - delta_rho[j]*dPot_dPos_trans[kj_idx]);
            total_forces[k] -= 0.5*(delta_rho[k]*delta_rho[k] + delta_rho[j]*delta_rho[j])*dDamp_dPos[kj_idx];
        }
        total_forces[k] += frozen_forces[k];
        //printf(" %.8f  %.8f  %.8f \n", dJ_dx[i*n_sites + i], dJ_dy[i*n_sites + i], dJ_dz[i*n_sites + i]);
    }
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
    printf(" Number of elements: %d \n", (int)params->size());
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
    return vec_d(J_mat);
}

vec_d FlucDens::get_rho_pot_vec()
{
    return vec_d(pot_vec);
}

vec_d FlucDens::get_delta_rho()
{
    return delta_rho;
}

double FlucDens::get_frozen_energy()
{
    return total_energies.frz;
}

double FlucDens::get_polarization_energy()
{
    return total_energies.pol;
}

FlucDensEnergies FlucDens::get_energies()
{
    return total_energies;
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