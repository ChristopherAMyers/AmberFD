#include "DispersionPauli.h"

using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::set;
using std::map;

DispersionPauli::DispersionPauli(const int num_sites, const int* nuclei_in, const double* exponents_in, const double* radii_in)
{
    n_sites = num_sites;

    //  copy over parameters
    nuclei.assign(nuclei_in, nuclei_in + n_sites);
    pauli_exponents.assign(exponents_in, exponents_in + n_sites);
    pauli_radii.assign(radii_in, radii_in + n_sites);

    //  convert radii over to pauli coefficients
    pauli_coeff.resize(pauli_radii.size(), 0.0);
    for (size_t i=0; i < pauli_radii.size(); i ++)
        pauli_coeff[i] = radii_to_coeff(pauli_radii[i], pauli_exponents[i]);
    // coeff12 = np.sqrt((4.184/2625.5009)*np.exp(exp_list12*radii*ang2bohr))

    //  initialize exclusions
    exclusions.resize(num_sites);

    //  default dispersion params
    disp_a1 = 0.294;
    disp_a2 = 2.266;
    disp_s6 = 1.200;
    
    //  default C6 coefficients
    C6_map = {
        {1, 4.30},
        {6, 6.55},
        {7, 6.17},
        {8, 5.62}
    };
    set_all_C6_coeff();

    //  default vdW radii
    vdw_radii_map = {
        {1, 1.001*ANG2BOHR},
        {6, 1.452*ANG2BOHR},
        {7, 1.397*ANG2BOHR},
        {8, 1.342*ANG2BOHR}
    };
    set_all_vdw_radii();

    self_forces.resize(n_sites, Vec3());
    secondary_radii_map = {
        {0, 2.1*0.35*ANG2BOHR},
        {1, 0.0*0.35*ANG2BOHR},
        {6, 2.8*0.35*ANG2BOHR},
        {7, 2.4*0.35*ANG2BOHR},
        {8, 2.2*0.35*ANG2BOHR}
    };
    // secondary_radii_map = {
    //     {0, 2.1*0.35*ANG2BOHR},
    //     {1, 0.0*0.35*ANG2BOHR},
    //     {6, 2.8*0.35*ANG2BOHR},
    //     {7, 2.4*0.35*ANG2BOHR},
    //     {8, 2.2*0.35*ANG2BOHR}
    // };
    secondary_exp = 8.0/ANG2BOHR;
    set_all_secondary_radii();
    use_secondary_radii = true;


    use_two_site_repulsion = false;
    two_site_indicies_set = false;

}

DispersionPauli::~DispersionPauli(){}

void DispersionPauli::set_all_secondary_radii()
{
    secondary_radii.resize(n_sites, 0.0);
    for(int i = 0; i < n_sites; i++)
        secondary_radii[i] = secondary_radii_map[nuclei[i]]; 
}

void DispersionPauli::set_all_vdw_radii()
{
    vdw_radii.resize(n_sites, 0.0);
    for(int i = 0; i < n_sites; i++)
        vdw_radii[i] = vdw_radii_map[nuclei[i]];
}

void DispersionPauli::set_all_C6_coeff()
{
    C6_coeff.resize(n_sites, 0.0);
    for(int i = 0; i < n_sites; i++)
        C6_coeff[i] = C6_map[nuclei[i]];
}


//****************************************************************
//                   Parameter Specification
//****************************************************************
void DispersionPauli::set_dispersion_params(double s6, double a1, double a2)
{
    disp_s6 = s6;
    disp_a1 = a1;
    disp_a2 = a2;
}

void DispersionPauli::get_dispersion_params(double &s6, double &a1, double &a2)
{
    s6 = disp_s6;
    a1 = disp_a1;
    a2 = disp_a2;
}

void DispersionPauli::set_vdw_radii(map_id nuclei2radiiMap)
{
    vdw_radii_map = nuclei2radiiMap;
    set_all_vdw_radii();
}

map_id DispersionPauli::get_vdw_radii_map()
{
    return vdw_radii_map;
}

vec_d DispersionPauli::get_vdw_radii()
{
    return vdw_radii;
}

void DispersionPauli::set_C6_map(map_id nucleiToC6Map){
    C6_map = nucleiToC6Map;
    set_all_C6_coeff();
}

map_id DispersionPauli::get_C6_map()
{
    return C6_map;
}

vec_d DispersionPauli::get_C6_coeff()
{
    return C6_coeff;
}

double DispersionPauli::radii_to_coeff(double radii, double exponent)
{
    if (radii <= 0)
        return 0.0;
    else
        return sqrt(kcal * exp(exponent*radii));
}

void DispersionPauli::set_pauli_radii(vec_d radii_list)
{
    if ((int)radii_list.size() != n_sites)
        throw std::runtime_error("radi_list length does not equal n_sites");
    pauli_radii.assign(radii_list.begin(), radii_list.end());
    for (size_t i=0; i < pauli_radii.size(); i ++)
        pauli_coeff[i] = radii_to_coeff(pauli_radii[i], pauli_exponents[i]);
}

void DispersionPauli::set_pauli_radii(int index, double radii)
{
    if (index > (n_sites - 1))
        throw std::out_of_range("radii index is greater than n_sites");
    if (index < 0)
        throw std::out_of_range("radii index must be greater than zero");
    pauli_radii[index] = radii;
    pauli_coeff[index] = radii_to_coeff(radii, pauli_exponents[index]);
}

vec_d DispersionPauli::get_pauli_radii()
{
    return pauli_radii;
}

void DispersionPauli::set_pauli_exp(vec_d exp_list)
{
    if ((int)exp_list.size() != n_sites)
        throw std::runtime_error("exp_list length does not equal n_sites");
    pauli_exponents.assign(exp_list.begin(), exp_list.end());
}

void DispersionPauli::set_pauli_exp(int index, double exponent)
{
    if (index > (n_sites - 1))
        throw std::out_of_range("exponent index is greater than n_sites");
    if (index < 0)
        throw std::out_of_range("exponent index must be greater than zero");
    pauli_exponents[index] = exponent;
}

vec_d DispersionPauli::get_pauli_exp()
{
    return pauli_exponents;
}
int DispersionPauli::get_num_sites()
{
    return n_sites;
}


//****************************************************************
//                      Exclusions
//****************************************************************
void DispersionPauli::create_exclusions_from_bonds(const vector<pair<int, int> > bonds, int bond_cutoff)
{
    vector<set<int> > exclusions = Nonbonded::calc_exclusions_from_bonds(bonds, bond_cutoff, n_sites);

    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i)
                add_exclusion(i, j);
}

void DispersionPauli::create_exclusions_from_fragment(const std::vector<int> frag_idx)
{
    for (auto idx_i: frag_idx)
        for(auto idx_j: frag_idx)
            add_exclusion(idx_i, idx_j);
}

void DispersionPauli::add_exclusion(const int index_i, const int index_j)
{
    if ((index_i > n_sites) || (index_j > n_sites) || (index_i < 0) || (index_j < 0))
    {   
        char buffer[100];
        sprintf(buffer, " Exception index pair (%d,%d) out of bounds", index_i, index_j);
        throw std::out_of_range(buffer);
    }
    exclusions[index_i].insert(index_j);
    exclusions[index_j].insert(index_i);
}

std::set<int> DispersionPauli::get_exclusions(const int particle1) const
{
    std::set<int> particles2;
    if (particle1 < 0 || particle1 >= (int) exclusions.size()) 
        throw "Index out of range";
    particles2 = exclusions[particle1];
    return particles2;
}


//****************************************************************
//                 Periodic Boundary Conditions
//****************************************************************
void DispersionPauli::set_use_PBC(const bool is_periodic)
{
    periodicity.is_periodic = is_periodic;
}
void DispersionPauli::set_use_PBC(const bool is_periodic, const double x, const double y, const double z)
{
    periodicity.set(is_periodic, x, y, z);
}

bool DispersionPauli::get_use_PBC()
{
    return periodicity.is_periodic;
}


//****************************************************************
//                          Energies
//****************************************************************
double DispersionPauli::get_pauli_energy()
{
    return total_pauli_energy;
}
double DispersionPauli::get_disp_energy()
{
    return total_disp_energy;
}

std::vector<vec_d> DispersionPauli::get_forces()
{
    std::vector<vec_d> rtn(n_sites, vec_d(3, 0.0));
    for(size_t i = 0; i < self_forces.size(); i++)
    {
        rtn[i][0] = self_forces[i][0];
        rtn[i][1] = self_forces[i][1];
        rtn[i][2] = self_forces[i][2];
    }
    return rtn;
}

void DispersionPauli::initialize()
{
    total_disp_energy = 0.0;
    total_pauli_energy = 0.0;
    std::fill(self_forces.begin(), self_forces.end(), Vec3(0, 0, 0));
}

double DispersionPauli::calc_energy(const vec_d &positions)
{
    size_t i, j;
    DeltaR deltaR;
    initialize();

    Energies energies;
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            if (periodicity.is_periodic)
                deltaR.getDeltaR(positions, (int)i*3, (int)j*3, periodicity);
            else
                deltaR.getDeltaR(positions, (int)i*3, (int)j*3);

            calc_one_pair(positions, deltaR, i, j, energies, self_forces);
            total_pauli_energy += energies.pauli;
            total_pauli_energy += energies.pauli_wall;
            total_disp_energy += energies.disp;
        }
    }

    return total_disp_energy + total_pauli_energy;
}

void DispersionPauli::add_force(vec_d &force, const Vec3 &dR)
{
    force[0] += dR[0];
    force[1] += dR[1];
    force[2] += dR[2];
}

Energies DispersionPauli::calc_one_pair(const vec_d &pos, int i, int j, std::vector<Vec3> &forces)
{
    Energies energies;
    DeltaR deltaR;
    if (periodicity.is_periodic)
        deltaR.getDeltaR(pos, (int)i*3, (int)j*3, periodicity);
    else
        deltaR.getDeltaR(pos, (int)i*3, (int)j*3);
    calc_one_pair(pos, deltaR, i, j, energies, forces);
    return energies;
}

double DispersionPauli::calc_one_pair(const vec_d &pos, DeltaR &deltaR, int i, int j, Energies& energies, std::vector<Vec3> &forces)
{
    energies.disp = 0.0;
    energies.pauli = 0.0;
    energies.pauli_wall = 0.0;
    if (deltaR.r > 30.0)
        return energies.disp + energies.pauli + energies.pauli_wall;
    if (std::find(exclusions[i].begin(), exclusions[i].end(), j) == exclusions[i].end())
    {      
        Vec3 dR_dPos = deltaR.dR*deltaR.r_inv;

        //  dispersion energy
        double radii = vdw_radii[i] + vdw_radii[j];
        double shift = disp_a1*radii + disp_a2;
        double shift2 = shift*shift;
        double shift6 = shift2*shift2*shift2;
        double r2 = deltaR.r2;
        double r6 = r2*r2*r2;
        double C6 = sqrt(C6_coeff[i]*C6_coeff[j]);
        double r6_shift_inv = 1/(r6 + shift6);

        energies.disp = -disp_s6*C6*r6_shift_inv;
        
        //  dispersion forces
        double dE_dR = -energies.disp * 6 * r6 * deltaR.r_inv * r6_shift_inv;
        forces[i] -=dE_dR*dR_dPos;
        forces[j] +=dE_dR*dR_dPos;

        //  Pauli energy
        double coeff = pauli_coeff[i]*pauli_coeff[j];
        double exponent = 0.5*(pauli_exponents[i] + pauli_exponents[j]);
        energies.pauli = coeff*exp(-exponent*deltaR.r);

        //  Pauli forces
        dE_dR = -energies.pauli*exponent;
        forces[i] -=dE_dR*dR_dPos;
        forces[j] +=dE_dR*dR_dPos;


        if (use_secondary_radii)
        {
            // backup secondary force energy
            double radii_2 = secondary_radii[i] + secondary_radii[j];
            double eng = kcal*exp(-secondary_exp*(deltaR.r - radii_2));
            energies.pauli_wall = eng;

            // backup secondary force force
            dE_dR = -eng*secondary_exp;
            forces[i] -=dE_dR*dR_dPos;
            forces[j] +=dE_dR*dR_dPos;
        }
    }
    return energies.disp + energies.pauli + energies.pauli_wall;
}

Vec3 DispersionPauli::get_perp_vector(const vec_d &pos, int i, int j, int k)
{
    // return vector normal to the plane formed by three sites
    Vec3 r1(pos[i*3], pos[i*3+1], pos[i*3+2]);
    Vec3 r2(pos[j*3], pos[j*3+1], pos[j*3+2]);
    Vec3 r3(pos[k*3], pos[k*3+1], pos[k*3+2]);

    Vec3 dx = r2 - r1;
    Vec3 dy = r3 - r1;
    Vec3 dz = dx.cross(dy);
    dz *= 1.0/sqrt(dz.dot(dz));
    // double normZ = sqrt(dz.dot(dz));
    // double invNorm = (normZ > 0.0 ? 1.0/normZ : 0.0);
    // dz *= invNorm;
    return dz;
}

void DispersionPauli::calc_two_site_repulsion(const vec_d &pos, DeltaR &deltaR, int i, int j, Energies& energies)
{
    double coeff = pauli_coeff[i]*pauli_coeff[j];
    double exponent = 0.5*(pauli_exponents[i] + pauli_exponents[j]);
    DeltaR dR;

    vector<Vec3> sites_i, sites_j;
    Vec3 r1(pos[i*3], pos[i*3+1], pos[i*3+2]);
    Vec3 r2(pos[j*3], pos[j*3+1], pos[j*3+2]);
    Vec3 dz;

    if (two_site_indicies[i].first != -1)
    {
        int idx_2 = two_site_indicies[i].first;
        int idx_3 = two_site_indicies[i].second;
        dz = get_perp_vector(pos, i, idx_2, idx_3);
        Vec3 r1_p = r1 + two_site_dist*dz;
        Vec3 r1_m = r1 - two_site_dist*dz;
        sites_i.push_back(r1_m);
        sites_i.push_back(r1_p);
    }
    else
        sites_i.push_back(r1);

    if (two_site_indicies[j].first != -1)
    {
        int idx_2 = two_site_indicies[j].first;
        int idx_3 = two_site_indicies[j].second;
        dz = get_perp_vector(pos, j, idx_2, idx_3);
        Vec3 r2_p = r2 + two_site_dist*dz;
        Vec3 r2_m = r2 - two_site_dist*dz;
        sites_j.push_back(r2_p);
        sites_j.push_back(r2_m);
    }
    else
        sites_j.push_back(r2);
    
    coeff /= ((double)sites_i.size()*(double)sites_j.size());
    
    // printf("NUM SITES: %d  %d \n", (int)sites_i.size(), (int)sites_j.size());
    for (auto &ri: sites_i)
    {
        for (auto &rj: sites_j)
        {
            if (periodicity.is_periodic)
                dR.getDeltaR(ri, rj, periodicity);
            else
                dR.getDeltaR(ri, rj);
            energies.pauli += coeff*exp(-exponent*dR.r);
            // printf("    DR:  %.5f \n", dR.r);
            // cout << "    " << ri << endl;
            // cout << "    " << rj << endl;
            // cout << "    " << dz << endl;
        }   
    }
    // printf("IN PAULI: %d  %d  %.5f \n", i, j, energies.pauli);
}

//****************************************************************
//                      Miscellaneous
//****************************************************************
void DispersionPauli::set_use_secondary_radii(bool use_radii)
{
    use_secondary_radii = use_radii;
}
void DispersionPauli::set_use_two_site_repulsion(bool on_off)
{
    use_two_site_repulsion = on_off;
}

void DispersionPauli::set_two_site_distance(double distance)
{
    two_site_dist = distance;
}

void DispersionPauli::create_repulsion_sites(double vertical_dist, const std::vector<std::pair<int, int>> &bonds)
{
    two_site_dist = vertical_dist;
    two_site_indicies.resize(n_sites, std::pair<int, int>(-1, -1));
    std::vector<vec_i> bonded_to;
    bonded_to.resize(n_sites);
    for(auto &pair: bonds)
    {
        bonded_to[pair.first].push_back(pair.second);
        bonded_to[pair.second].push_back(pair.first);
    }

    for(int i = 0; i < n_sites; i++)
    {
        //  exclude hydrogens
        if (nuclei[i] > 2)
        {
            vec_i bonds = bonded_to[i];
            
            //  only atoms with one or two bonds can have out-of-plane pauli sites
            if (bonds.size() == 2)
                two_site_indicies[i] = std::pair<int, int>(bonds[0], bonds[1]);
            else if (bonds.size() == 1)
            {
                int site_1 = bonds[0];
                int site_2 = bonded_to[site_1][0];
                if (site_2 == i)
                    site_2 = bonded_to[site_1][1];
                two_site_indicies[i] = std::pair<int, int>(site_1, site_2);
            }
        }
    }
}

