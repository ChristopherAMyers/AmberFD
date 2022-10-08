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
    remove_core = true;
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
    {
        if ((nuclei[i] > 2) && remove_core)
            nuclei[i] -= 2;
        frozen_pop.push_back(nuclei[i] - frozen_chg[i]);
    }

    //  param_data is for easy access to parameters my names
    param_data["frz_chg"] =    &frozen_chg;
    param_data["frz_pop"] =    &frozen_pop;
    param_data["frz_exp"] =    &frozen_exp;
    param_data["dyn_exp"] =    &dynamic_exp;
    param_data["nuclei"] =     &nuclei;

    //  for dynamic densities
    J_mat.resize(n_sites*n_sites, 0.0);
    dJ_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dPot_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dPot_dPos_trans.resize(n_sites*n_sites, Vec3(0, 0, 0));
    dDamp_dPos.resize(n_sites*n_sites, Vec3(0, 0, 0));
    pot_vec.resize(n_sites, 0.0);
    potential_mat.resize(n_sites*n_sites, 0.0);
    for (int i = 0; i < n_sites; i++)
        J_mat[i*n_sites + i] = dynamic_exp[i]*5/16;
    hardness.resize(n_sites, 0.0);

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
    dampening_type = DampType::Linear;
    thread_dampening.resize(Nonbonded::num_threads);
    for (int i = 0; i < Nonbonded::num_threads; i++)
        thread_dampening[i].resize(n_sites, 0.0);
    ct_coeff = 0.0;
    calc_forces = true;

    //  initialize energies
    initialize_calculation();

    //  cutoff info;
    use_cutoff = false;
    use_SR_cutoff = true;
    cutoff_distance = 0;
    SR_cutoff_coeff = SR_cutoff_a*log(SR_cutoff_pct_error) + SR_cutoff_b;

    //  external fields
    has_ext_field = false;
    ext_field_potential.resize(n_sites, 0.0);

    //  default minimizer
    solver = Solver::Global;
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
//****************************************************************
//                   individual site parameters 
//****************************************************************
void FlucDens::set_site_params(const int index, const double frz_chg_new, const double frz_exp_new, const double dyn_exp_new)
{
    frozen_chg[index] = frz_chg_new;
    frozen_exp[index] = frz_exp_new;
    dynamic_exp[index] = dyn_exp_new;
    frozen_pop[index] = nuclei[index] - frozen_chg[index];
}

void FlucDens::get_site_params(const int index, double &frz_chg, double &frz_exp, double &dyn_exp)
{
    frz_chg = frozen_chg[index];
    frz_exp = frozen_exp[index];
    dyn_exp = dynamic_exp[index];
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

void FlucDens::set_ct_coeff(const double coeff)
{
    ct_coeff = coeff;
}

double FlucDens::get_ct_coeff()
{
    return ct_coeff;
}

void FlucDens::set_dampening(double coeff, double exponent, DampType damp)
{
    damp_coeff = coeff;
    damp_exponent = exponent;
    dampening_type = damp;
}

void FlucDens::get_dampening(double &coeff, double &exponent)
{
    coeff = damp_coeff;
    exponent = damp_exponent;
}

void FlucDens::set_additional_hardness(const int index, const double value)
{
    hardness[index] = value;
}

void FlucDens::set_additional_hardness(vec_d values)
{
    if ((int)values.size() != n_sites)
        throw std::runtime_error("hardness array length does not equal n_sites");
    hardness.assign(values.begin(), values.end());
}

//****************************************************************
//**********  delta_rho - frozen exclusions and fragments ********
//****************************************************************
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

std::vector<vec_i> FlucDens::get_fragments()
{
    if (site_frag_ids_partitioned.size() != n_fragments)
    {
        site_frag_ids_partitioned.resize(n_fragments);
        for (int i = 0; i < n_sites; i++)
        {
            site_frag_ids_partitioned[site_frag_ids[i]].push_back(i);
        }
    }
    return site_frag_ids_partitioned;
}

int FlucDens::get_num_fragments()
{
    return n_fragments;
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

void FlucDens::set_frag_constraints(const bool constr_frags)
{
    /* choose to constrain total charge per fragment (True) or as a whole system (False) */
    use_frag_constraints = constr_frags;
}

int FlucDens::get_num_constraints()
{
    return constraints.size();
}

std::vector<std::vector<int>> FlucDens::get_constraints()
{
    return constraints;
}

void FlucDens::assign_constraints()
{
    if (use_frag_constraints)
    {
        if (n_fragments == 0)
            throw std::runtime_error("No fragments defined when fragment constraints are enabled");
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

void FlucDens::create_del_exclusions_from_fragment(const std::vector<int> frag_idx)
{
    for (auto del_i: frag_idx)
        for(auto frz_j: frag_idx)
            add_del_frz_exclusion(del_i, frz_j);
}


//****************************************************************
//                frozen - frozen exclusions 
//****************************************************************
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

std::set<int> FlucDens::get_frz_frz_exclusions(const int particle1) const
{
    if (particle1 < 0 || particle1 >= (int) exclusions_frz_frz.size()) 
        throw "Index out of range";
    // std::vector<int> particles2;
    // particles2.resize(exclusions_frz_frz[particle1].size());
    // std::copy(exclusions_frz_frz[particle1].begin(), exclusions_frz_frz[particle1].end(), particles2.begin());
    // int i = 0;
    // for(auto &x: exclusions_frz_frz[particle1])
    // {
    //     printf("%d  %d  %d \n", (int)i, particles2[i], x);
    //     i++;
    // }
    // return particles2;
    return exclusions_frz_frz[particle1];
}

int FlucDens::get_num_frz_frz_exclusions() const
{
    return exclusions_frz_frz.size(); 
}

void FlucDens::create_frz_exclusions_from_bonds(const vector<pair<int, int> > bonds, int bond_cutoff)
{
    vector<std::set<int> > exclusions = Nonbonded::calc_exclusions_from_bonds(bonds, bond_cutoff, n_sites);
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i)
                add_frz_frz_exclusion(i, j);
}

//****************************************************************
//                  External fields and dipoles 
//****************************************************************
void FlucDens::set_external_field(double field_x, double field_y, double field_z)
{
    has_ext_field = true;
    ext_field[0] = field_x;
    ext_field[1] = field_y;
    ext_field[2] = field_z;
}

void FlucDens::apply_field_to_system(const vec_d &coords)
{
    // if (periodicity.is_periodic)
    //     throw std::runtime_error("External electric field is not (yet) vallid with periodic systems");
    if (has_ext_field)
    {
        // printf("Applying external field %.5f  %.5f  %.5f \n", ext_field[0], ext_field[1], ext_field[2]);
        for(int i = 0; i < n_sites; i++)
        {
            ext_field_potential[i] = -ext_field[0]*coords[i*3] - ext_field[1]*coords[i*3 + 1] - ext_field[2]*coords[i*3 +2];
        }
    }
}

std::vector<Vec3> FlucDens::get_dipoles(const vec_d &coords)
{
    /* return the dipole moment from the nuclei, frozen electrons, 
    and dynamic electrons */
    std::vector<Vec3> dipoles(4, Vec3(0.0, 0.0, 0.0));
    for (int i = 0; i < n_sites; i++)
    {
        Vec3 pos(coords[i*3 + 0], coords[i*3 + 1], coords[i*3 + 2]);
        dipoles[DensityType::Nuclei] += pos*nuclei[i];
        dipoles[DensityType::Frozen] -= pos*frozen_pop[i];
        dipoles[DensityType::Delta] -= pos*delta_rho[i];
    }
    dipoles[DensityType::All] = dipoles[DensityType::Nuclei] + dipoles[DensityType::Frozen] + dipoles[DensityType::Delta];
    return dipoles;
}

vec_d FlucDens::get_dipole(const vec_d &coords, DensityType density_dype=DensityType::All)
{
    std::vector<Vec3> dipoles = get_dipoles(coords);
    vec_d result = {dipoles[density_dype][0], dipoles[density_dype][1], dipoles[density_dype][2]};
    return result;
}

double FlucDens::calc_frz_ext_field_energy(const vec_d &positions, std::vector<Vec3> &forces)
{
    double energy = 0.0;
    if (has_ext_field)
    {
        for(int i = 0; i < n_sites; i ++)
            forces[i] += ext_field*frozen_chg[i];

        std::vector<Vec3> dipoles = get_dipoles(positions);
        Vec3 dipole_total =  dipoles[DensityType::Nuclei] + dipoles[DensityType::Frozen];
        energy = -dipole_total.dot(ext_field);
    }
    return energy;
}

//****************************************************************
//                  Periodic boundary conditions
//****************************************************************
void FlucDens::set_use_PBC(const bool is_periodic)
{
    periodicity.is_periodic = is_periodic;
}

void FlucDens::set_use_PBC(const bool is_periodic, const double x, const double y, const double z)
{
    periodicity.set(is_periodic, x, y, z);
}

bool FlucDens::get_use_PBC()
{
    return periodicity.is_periodic;
}

//****************************************************************
//                      Cutoff distances
//****************************************************************
void FlucDens::set_use_cutoff(bool useCutoff)
{
    use_cutoff = useCutoff;
}

bool FlucDens::get_use_cutoff()
{
    return use_cutoff;
}
void FlucDens::set_cutoff_distance(double distance_in_nm)
{
    cutoff_distance = distance_in_nm*10*ANG2BOHR;
}

double FlucDens::get_cutoff_distance()
{
    return cutoff_distance/(10*ANG2BOHR);
}

void FlucDens::set_use_SR_cutoff(bool yes_no)
{
    use_SR_cutoff = yes_no;
}

bool FlucDens::get_use_SR_cutoff()
{
    return use_SR_cutoff;
}

bool FlucDens::use_SR_approx(double r, double a, double b)
{
    if (!use_SR_cutoff)
        return false;
    double alpha = std::max(a, b);
    double max_dist = SR_cutoff_coeff/alpha;
    return (r > max_dist);
}

//****************************************************************
//                      Energies and Forces
//****************************************************************
double FlucDens::get_frozen_energy()
{
    return total_energies.frz;
}

double FlucDens::get_polarization_energy()
{
    return total_energies.pol;
}

double FlucDens::get_ct_energy()
{
    return total_energies.vct;
}

FlucDensEnergies FlucDens::get_energies()
{
    return total_energies;
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

double FlucDens::calc_overlap(const vec_d &coords, DensityType density_type, const vec_i &indices_in)
{
    size_t i, j;
    double overlap = 0.0;
    double r, r2, inv_r;
    double pop_a, pop_b, a, b, exp_ar, exp_br;

    vec_i indices(indices_in);
    
    if (indices.size() == 0)
    {
        indices.resize(n_sites);
        for(int i = 0; i < n_sites; i++)
        {
            indices[i] = i;
        }
    }

    vec_d *pops, *exponents;
    if(density_type == DensityType::Delta)
    {
        pops = &delta_rho;
        exponents = &dynamic_exp;
    }
    else if(density_type == DensityType::Frozen)
    {
        pops = &frozen_pop;
        exponents = &frozen_exp;
    }
    else if(density_type == DensityType::All)
    {
        throw std::invalid_argument("Overlap for ALL not implimented yet");
    }

    for(int i: indices)
    {
        a = (*exponents)[i];
        pop_a = (*pops)[i];
        

        for(int j: indices)
        {
            if (j < i)
                continue;

            r2 = dot3Vec(coords, (int)j*3, (int)i*3);
            r = sqrt(r2);
            inv_r = 1/r;
            b = (*exponents)[j];
            pop_b = (*pops)[j];
            
            //if (std::find(exclusions_frz_frz[i].begin(), exclusions_frz_frz[i].end(), j) == exclusions_frz_frz[i].end())
            {
                if (use_SR_approx(r, a, b))
                { }
                else if (true)
                {
                    //  compute exponentials only once
                    exp_br = exp(-b*r);
                    exp_ar = exp(-a*r);
                    double S_12 = frz_frz_overlap(inv_r, a, b, exp_ar, exp_br);
                    overlap += S_12*pop_a*pop_b;
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
        //double r = 1/inv_r;
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
    //fill(damp_sum.begin(), damp_sum.end(), 0.0);
    fill(pot_vec.begin(), pot_vec.end(), 0.0);
    fill(J_mat.begin(), J_mat.end(), 0.0);
    fill(dJ_dPos.begin(), dJ_dPos.end(), Vec3(0, 0, 0));
    fill(dPot_dPos.begin(), dPot_dPos.end(), Vec3(0, 0, 0));
    fill(dPot_dPos_trans.begin(), dPot_dPos_trans.end(), Vec3(0, 0, 0));
    fill(dDamp_dPos.begin(), dDamp_dPos.end(), Vec3(0, 0, 0));
    fill(frozen_forces.begin(), frozen_forces.end(), Vec3(0, 0, 0));
    fill(total_forces.begin(), total_forces.end(), Vec3(0, 0, 0));
    total_energies.reset();
    if (solver == Solver::Global)
        fill(delta_rho.begin(), delta_rho.end(), 0.0);

    if ((int)cutoff_distance == 0 && periodicity.is_periodic)
    {
        double min_dist = std::min(periodicity.box_size[0], periodicity.box_size[1]);
        min_dist = std::min(min_dist, periodicity.box_size[2]);
        cutoff_distance = 0.5*0.99999*min_dist;
        cutoff_distance = std::min(40.0, cutoff_distance);
    }
    if (Nonbonded::num_threads != (int)thread_dampening.size())
    {
        thread_dampening.resize(Nonbonded::num_threads);
        for (int i = 0; i < Nonbonded::num_threads; i++)
            thread_dampening[i].resize(n_sites);
    }
    for (int i = 0; i < Nonbonded::num_threads; i++)
        fill(thread_dampening[i].begin(), thread_dampening[i].end(), 0.0);
    

}

double FlucDens::calc_energy(const vec_d &positions, bool calc_frz, bool calc_pol)
{
    initialize_calculation();
    if (!calc_pol && !calc_frz)
        return 0.0;

    size_t i, j;
    Energies energies;
    DeltaR deltaR;
    energies.zero();
    double total = 0.0;

    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances and inverse distances
            if (periodicity.is_periodic)
                deltaR.getDeltaR(positions, (int)i*3, (int)j*3, periodicity);
            else
                deltaR.getDeltaR(positions, (int)i*3, (int)j*3);
            
            calc_one_electro(deltaR, i, j, calc_pol, calc_frz, energies, total_forces);
            total_energies.frz += energies.frz;
            total_energies.nuc_nuc += energies.nuc_nuc;
            total_energies.elec_elec += energies.elec_elec;
            total_energies.elec_nuc += energies.elec_nuc;
        }
    }
    if (calc_pol)
    {
        solve_minimization(total_forces);
    }
    total += calc_frz_ext_field_energy(positions, total_forces);
    // if (has_ext_field)
    // {
    //     for(int i = 0; i < n_sites; i ++)
    //         total_forces[i] += ext_field*(frozen_chg[i] - delta_rho[i]);
    //     std::vector<Vec3> dipoles = get_dipoles(positions);
    //     Vec3 dipole_total =  dipoles[0] + dipoles[1] + dipoles[2];
    //     total  -= dipole_total.dot(ext_field);
    // }
    total += total_energies.total();

    return total;
}

Energies FlucDens::calc_one_frozen(const vec_d &positions, int i, int j)
{
    Energies eng_out;
    DeltaR deltaR;
    if (periodicity.is_periodic)
        deltaR.getDeltaR(positions, (int)i*3, (int)j*3, periodicity);
    else
        deltaR.getDeltaR(positions, (int)i*3, (int)j*3);
    calc_one_electro(deltaR, i, j, false, true, eng_out, total_forces);
    return eng_out;
}

void FlucDens::calc_one_electro(DeltaR &deltaR, int i, int j, bool calc_pol, bool calc_frz, Energies& energies, std::vector<Vec3> &forces, int thread_num)
{
    /*  calculate a single electrostatic interaction between a pair of atoms */
    
    double frozen_energy = 0.0;
    double E_ee = 0.0, E_eZ1 = 0.0, E_eZ2 = 0.0, E_ZZ = 0.0;

    bool within_cutoff = !use_cutoff || deltaR.r < cutoff_distance;

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

    if (use_SR_approx(r, a_frz, b_frz)
    && use_SR_approx(r, a_del, b_del) )
    {            
        Vec3 dR = deltaR.dR*inv_r;
        double dE_dR = -inv_r*inv_r;

        if (calc_pol and within_cutoff)
        {
            //  delta_rho - delta_rho
            J_mat[idx_ij] = inv_r;
            J_mat[idx_ji] = inv_r;
            dJ_dPos[idx_ij] =  dE_dR * dR;
            dJ_dPos[idx_ji] = -dE_dR * dR;

            //  delta_rho - frozen_rho
            if (
                    (exclusions_del_frz[i].find(j) == exclusions_del_frz[i].end())
                    && within_cutoff
                )
            {
                potential_mat[idx_ij] = -frozen_chg[j]*inv_r;
                potential_mat[idx_ji] = -frozen_chg[i]*inv_r;

                double dP_dR_ij = -frozen_chg[j]*dE_dR;
                double dP_dR_ji = -frozen_chg[i]*dE_dR;
                dPot_dPos[idx_ij] =  dP_dR_ij*dR;
                dPot_dPos[idx_ji] = -dP_dR_ji*dR;
                dPot_dPos_trans[idx_ji] =  dP_dR_ij*dR;
                dPot_dPos_trans[idx_ij] = -dP_dR_ji*dR;

                //  polarization dampening
                // double dampening = exp(-damp_exponent*r)*0;
                // thread_dampening[thread_num][i] += dampening*damp_coeff;
                // thread_dampening[thread_num][j] += dampening*damp_coeff;
                // double dDamp_dR = -damp_exponent*damp_coeff*dampening;
                // dDamp_dPos[idx_ij] =  dDamp_dR*dR;
                // dDamp_dPos[idx_ji] = -dDamp_dR*dR;
            }
        }

        //  frozen_rho - frozen_rho
        if (calc_frz && within_cutoff)
        if (exclusions_frz_frz[i].find(j) == exclusions_frz_frz[i].end())
        //if (std::find(exclusions_frz_frz[i].begin(), exclusions_frz_frz[i].end(), j) == exclusions_frz_frz[i].end())
        {
            E_ZZ = nuclei[i]*nuclei[j]*inv_r;
            E_ee = frozen_pop[i]*frozen_pop[j]*inv_r;
            E_eZ1 = -frozen_pop[i]*nuclei[j]*inv_r;
            E_eZ2 = -frozen_pop[j]*nuclei[i]*inv_r;

            Vec3 dR_norm = deltaR.dR*inv_r;
            Vec3 d_ee_dR = frozen_pop[i]*frozen_pop[j]*dE_dR*dR_norm;
            Vec3 d_eZ_dR1 = -frozen_pop[i]*nuclei[j]*dE_dR*dR_norm;
            Vec3 d_eZ_dR2 = -frozen_pop[j]*nuclei[i]*dE_dR*dR_norm;
            Vec3 d_ZZ_dR =  nuclei[i]*nuclei[j]*dE_dR*dR_norm;

            Vec3 force_i = -(d_ee_dR + d_eZ_dR1 + d_eZ_dR2 + d_ZZ_dR);
            forces[i] += force_i;
            forces[j] -= force_i;

            frozen_energy =  (E_ee + E_eZ1) + (E_ZZ + E_eZ2);
        }

        // skip dampening
    }
    else
    {
        const double exp_ar_del = exp(-a_del*r);
        const double exp_ar_frz = exp(-a_frz*r);
        const double exp_br_frz = exp(-b_frz*r);
        const double exp_br_del = exp(-b_del*r);

        if (calc_pol and within_cutoff)
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
            if (
                (exclusions_del_frz[i].find(j) == exclusions_del_frz[i].end())
                //(std::find(exclusions_del_frz[i].begin(), exclusions_del_frz[i].end(), j) == exclusions_del_frz[i].end())
                && within_cutoff
            )
            {
                double del_frz_ee_ddR, frz_del_ee_ddR, del_a_nuc_b_ddR, del_b_nuc_a_ddR;
                const double del_frz_ee  = elec_elec_energy(inv_r, a_del, b_frz, exp_ar_del, exp_br_frz, del_frz_ee_ddR);
                const double frz_del_ee  = elec_elec_energy(inv_r, a_frz, b_del, exp_ar_frz, exp_br_del, frz_del_ee_ddR);
                const double del_a_nuc_b = elec_nuclei_energy(inv_r, a_del, exp_ar_del, del_a_nuc_b_ddR);
                const double del_b_nuc_a = elec_nuclei_energy(inv_r, b_del, exp_br_del, del_b_nuc_a_ddR);
                
                potential_mat[idx_ij] = frozen_pop[j]*del_frz_ee + nuclei[j]*del_a_nuc_b;
                potential_mat[idx_ji] = frozen_pop[i]*frz_del_ee + nuclei[i]*del_b_nuc_a;
                
                double dP_dR_ij = frozen_pop[j]*del_frz_ee_ddR + nuclei[j]*del_a_nuc_b_ddR;
                double dP_dR_ji = frozen_pop[i]*frz_del_ee_ddR + nuclei[i]*del_b_nuc_a_ddR;
                dPot_dPos[idx_ij] = dP_dR_ij*dR;
                dPot_dPos[idx_ji] = -dP_dR_ji*dR;
                dPot_dPos_trans[idx_ji] = dP_dR_ij*dR;
                dPot_dPos_trans[idx_ij] = -dP_dR_ji*dR;

                //  polarization dampening
                
                //  proportional to frozen density
                // damp_coeff = 0.0;
                double damp_i = pow(frozen_pop[j]*b_frz*b_frz*b_frz/(8*M_PI)*exp_br_frz, damp_exponent)*damp_coeff;
                double damp_j = pow(frozen_pop[i]*a_frz*a_frz*a_frz/(8*M_PI)*exp_ar_frz, damp_exponent)*damp_coeff;
                thread_dampening[thread_num][i] += damp_i;
                thread_dampening[thread_num][j] += damp_j;

                dDamp_dPos[idx_ij] = -damp_i*damp_exponent*b_frz*dR;
                dDamp_dPos[idx_ji] =  damp_j*damp_exponent*a_frz*dR;


                // double dampening = exp(-damp_exponent*r);
                // thread_dampening[thread_num][i] += dampening*damp_coeff*(frozen_pop[j] > 0);
                // thread_dampening[thread_num][j] += dampening*damp_coeff*(frozen_pop[i] > 0);
                // thread_dampening[thread_num][i] += dampening*damp_coeff;
                // thread_dampening[thread_num][j] += dampening*damp_coeff;
                // double dDamp_dR = -damp_exponent*damp_coeff*dampening;
                // dDamp_dPos[idx_ij] =  dDamp_dR*dR;
                // dDamp_dPos[idx_ji] = -dDamp_dR*dR;



            }
        }
        //  frozen_rho - frozen_rho
        if (calc_frz && within_cutoff)
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
            forces[i] += force_i;
            forces[j] -= force_i;

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

void FlucDens::solve_minimization(std::vector<Vec3> &forces)
{
    /*  solve the block matrix equation Ax = b with
        |J   , C  || Δρ |     | -ϕ |
        |C^T , 0  || λ  |  =  | -D |
        A is composed of the Coulomb and constraint matrices.
        x is the solution for the change in electron density and Lagrange multipliers
        b is the vector of the negative potential and resulting constraints
    */
    assign_constraints();
    int n_constr = (int)constraints.size();
    int dim = n_sites + n_constr;
    int info;
    vec_d A_mat(dim*dim, 0.0);
    vec_d B_vec(dim, 0);

    //   //  copy over dampening
    // for (int i = 0; i < n_sites; i++)
    //     J_mat[i*n_sites + i] = dynamic_exp[i]*5/16 + damp_sum[i]*damp_coeff;
    
    for (int i = 0; i < n_sites; i++)
        J_mat[i*n_sites + i] = dynamic_exp[i]*5/16 + hardness[i];

    for (int j = 0; j < n_sites; j++)
        pot_vec[j] = 0.0;
    
    //  copy over dampening and zero out potentials
    for (int i = 0; i < Nonbonded::num_threads; i++)
    {
        for (int j = 0; j < n_sites; j++)
        {
            if (dampening_type == Quadratic)
                J_mat[j*n_sites + j] += thread_dampening[i][j];
            else
                pot_vec[j] += thread_dampening[i][j];
        }
    }

    vec_d neg_pot(n_sites, 0.0);
    for (int i = 0; i < n_sites; i++)
    {
        //pot_vec[i] = 0.0;
        for(int j = 0; j < n_sites; j++)
        {
            pot_vec[i] += potential_mat[i*n_sites + j];
        }
        if (has_ext_field)
            pot_vec[i] -= ext_field_potential[i]/(1+ct_coeff);

        neg_pot[i] = -pot_vec[i];        
    }

    //vec_d neg_pot = pot_vec;
    //cblas_dscal(n_sites, -1, &neg_pot[0], 1);   //negate potential (right hand side of matrix equation)
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

    if (solver == Solver::Global)
    {
        A_mat_save = A_mat;
        B_vec_save = B_vec;
        
        //  solve matrix equation A_mat * x = B_vec for x
        int ipiv[dim];
        //openblas_set_num_threads(1);
        info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim, 1, &A_mat[0], dim, ipiv, &B_vec[0], dim);
        delta_rho.assign(B_vec.begin(), B_vec.begin() + n_sites);
    }
    else{
        std::vector<vec_i> fragments = get_fragments();
        divide_and_conquer.solve(J_mat, pot_vec, fragments, delta_rho, dynamic_exp);
    }

    //  simplified version of polarization energy
    if (false)
    {
        //  Since all constraints are linear in delta_rho and result in sums equal to zero,
        //  the polarization energy is simply 0.5 * delta_rho^T * pot_vec
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

    //#pragma omp parallel for
    for(int k = 0; k < n_sites; k++)
    {
        for(int j = 0; j < n_sites; j++)
        {
            int kj_idx = k*n_sites + j;
            //  force due to delta-delta interactions
            forces[k] -= delta_rho[k]*dJ_dPos[kj_idx]*delta_rho[j]*(1+ct_coeff);
            //  force due to frozen-delta interactions
            forces[k] -= (delta_rho[k]*dPot_dPos[kj_idx] - delta_rho[j]*dPot_dPos_trans[kj_idx])*(1+ct_coeff);

            //  force due to polarization dampening
            if (dampening_type == Quadratic)
                forces[k] -= 0.5*(delta_rho[k]*delta_rho[k] + delta_rho[j]*delta_rho[j])*dDamp_dPos[kj_idx]*(1+ct_coeff);
            else
                // TODO: encorporate a dDamp_dPos_transpose array for quicker access
                forces[k] -= (delta_rho[k]*dDamp_dPos[kj_idx] - delta_rho[j]*dDamp_dPos[j*n_sites + k])*(1+ct_coeff);
        }
    }

    total_energies.vct = ct_coeff*total_energies.pol;
}

std::vector<std::string> FlucDens::get_param_names()
{
    std::vector<std::string> keys;
    for(std::map<std::string, vec_d* >::iterator it = param_data.begin(); it != param_data.end(); ++it)
    {
        keys.push_back(it->first);
    }
    return keys;
}

vec_d FlucDens::get_params_by_name(const std::string param_name)
{
    if (param_data.count(param_name) > 0)
        return vec_d(param_data[param_name]->begin(), param_data[param_name]->end());
    else
        return vec_d(0);
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

vec_d FlucDens::calc_density(const vec_d &points, const vec_d &pos, DensityType density_type)
{
    vec_d density(int(points.size()/3), 0.0);
    vec_d frz_coeff(frozen_exp);
    vec_d dyn_coeff(frozen_exp);
    for(size_t i = 0; i < n_sites; i++)
    {
        frz_coeff[i] = std::pow(frozen_exp[i], 3)/(8*M_PI)*frozen_pop[i];
        dyn_coeff[i] = std::pow(dynamic_exp[i], 3)/(8*M_PI)*delta_rho[i];
        //printf("DYN COEFF %10.5e\n", dyn_coeff[i]);
    }
    for(size_t i = 0; i < int(points.size()/3); i++)
    {
        for(size_t j = 0; j < n_sites; j++)
        {
            double dx = points[i*3 + 0] - pos[j*3 + 0];
            double dy = points[i*3 + 1] - pos[j*3 + 1];
            double dz = points[i*3 + 2] - pos[j*3 + 2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double dr = sqrt(dr2);

            if(density_type == DensityType::All)
            {
                double frz = frz_coeff[j]*std::exp(-frozen_exp[j]*dr);
                double dyn = dyn_coeff[j]*std::exp(-dynamic_exp[j]*dr);
                density[i] += frz + dyn;
            }
            else if(density_type == DensityType::Frozen)
            {
                density[i] += frz_coeff[j]*std::exp(-frozen_exp[j]*dr);
            }
            else if(density_type == DensityType::Delta)
            {
                density[i] += dyn_coeff[j]*std::exp(-dynamic_exp[j]*dr);
            }
        }

    }

    return density;
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

void FlucDens::set_calc_forces(bool calculate_forces)
{
    calc_forces = calculate_forces;
}

double FlucDens::get_total_time()
{
    return total_time;
}
void FlucDens::set_solver(Solver solver_in)
{
    solver = solver_in;
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

