#include "AmberFD.h"

AmberFD::AmberFD()
{}

AmberFD::AmberFD(const int n_sites)
{
    nuclei.reserve(n_sites);
    frz_chg.reserve(n_sites);
    frz_exp.reserve(n_sites);
    dyn_exp.reserve(n_sites);
    pauli_coeff.reserve(n_sites);
    pauli_exp.reserve(n_sites);
}
AmberFD::~AmberFD()
{
    //delete[] flucDens;
    //delete[] dispersionPauli;
}

void AmberFD::add_particle(vec_d parameters)
{
    if (parameters.size() != 6)
        throw std::invalid_argument(" AmberFD particle must use 6 parameters");

    
}

void AmberFD::add_particle(ParticleInfo parameters)
{
    nuclei.push_back(parameters.nuclei);
    frz_exp.push_back(parameters.frz_exp);
    frz_chg.push_back(parameters.frz_chg);
    dyn_exp.push_back(parameters.dyn_exp);
    pauli_exp.push_back(parameters.pauli_exp);
    pauli_coeff.push_back(parameters.pauli_coeff);
}

void AmberFD::create_fluc_dens_force()
{
    vec_d nuc(nuclei.begin(), nuclei.end());
    flucDens = std::shared_ptr<FlucDens>(new FlucDens(n_sites, &frz_chg[0], &nuc[0], &frz_exp[0], &dyn_exp[0]));
}

void AmberFD::create_disp_pauli_force()
{

}