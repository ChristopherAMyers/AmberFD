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
}

DispersionPauli::~DispersionPauli(){}

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

void DispersionPauli::set_dispersion_params(double s6, double a1, double a2)
{
    disp_s6 = s6;
    disp_a1 = a1;
    disp_a2 = a2;
}
void DispersionPauli::set_vdw_radii(map_id nuclei2radiiMap)
{
    vdw_radii_map = nuclei2radiiMap;
    set_all_vdw_radii();
}
void DispersionPauli::set_C6_map(map_id nucleiToC6Map){
    C6_map = nucleiToC6Map;
    set_all_C6_coeff();
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
        throw std::out_of_range("coeff index is greater than n_sites");
    if (index < 0)
        throw std::out_of_range("coeff index must be greater than zero");
    pauli_radii[index] = radii;
    pauli_coeff[index] = radii_to_coeff(radii, pauli_exponents[index]);
}
void DispersionPauli::set_pauli_exp(vec_d exp_list)
{
    if ((int)exp_list.size() != n_sites)
        throw std::runtime_error("coeff_list length does not equal n_sites");
    pauli_exponents.assign(exp_list.begin(), exp_list.end());
}
void DispersionPauli::set_pauli_exp(int index, double exponent)
{
    if (index > (n_sites - 1))
        throw std::out_of_range("coeff index is greater than n_sites");
    if (index < 0)
        throw std::out_of_range("coeff index must be greater than zero");
    pauli_exponents[index] = exponent;
}
void DispersionPauli::get_dispersion_params(double &s6, double &a1, double &a2)
{
    s6 = disp_s6;
    a1 = disp_a1;
    a2 = disp_a2;
}

map_id DispersionPauli::get_vdw_radii()
{
    return vdw_radii_map;
}
map_id DispersionPauli::get_C6_map()
{
    return C6_map;
}
vec_d DispersionPauli::get_pauli_radii()
{
    return pauli_radii;
}
vec_d DispersionPauli::get_pauli_exp()
{
    return pauli_exponents;
}
vec_d DispersionPauli::get_C6_coeff()
{
    return C6_coeff;
}
double DispersionPauli::get_pauli_energy()
{
    return total_pauli_energy;
}
double DispersionPauli::get_disp_energy()
{
    return total_disp_energy;
}
int DispersionPauli::get_num_sites()
{
    return n_sites;
}

double DispersionPauli::calc_energy(const vec_d &positions)
{
    size_t i, j;
    total_disp_energy = 0.0;
    total_pauli_energy = 0.0;

    Energies energies;
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances data
            double deltaR[Nonbonded::RMaxIdx];
            Nonbonded::calc_dR(positions, (int)j*3, (int)i*3, deltaR);

            calc_one_pair(deltaR, i, j, energies);
            total_pauli_energy += energies.pauli;
            total_disp_energy += energies.disp;
        }
    }

    return total_disp_energy + total_pauli_energy;
}

double DispersionPauli::calc_one_pair(double *deltaR, int i, int j, Energies& energies)
{
    double pair_disp = 0.0;
    double pair_pauli = 0.0;
    if (std::find(exclusions[i].begin(), exclusions[i].end(), j) == exclusions[i].end())
    {      
        //  dispersion energy
        double radii = vdw_radii[i] + vdw_radii[j];
        double shift = disp_a1*radii + disp_a2;
        double shift2 = shift*shift;
        double shift6 = shift2*shift2*shift2;
        double r2 = deltaR[Nonbonded::R2Idx];
        double r6 = r2*r2*r2;
        double C6 = sqrt(C6_coeff[i]*C6_coeff[j]);
        pair_disp = -disp_s6*C6/(r6 + shift6);

        //  Pauli energy
        double coeff = pauli_coeff[i]*pauli_coeff[j];
        double exponent = 0.5*(pauli_exponents[i] + pauli_exponents[j]);
        pair_pauli = coeff*exp(-exponent*deltaR[Nonbonded::RIdx]);

        //printf("NEW: %d, %d, %.10f\n", i, j, pair_pauli);
        //  return and update totals
        // total_disp_energy += pair_disp;
        // total_pauli_energy += pair_pauli;
    }
    energies.disp = pair_disp;
    energies.pauli = pair_pauli;
    return pair_disp + pair_pauli;
}

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