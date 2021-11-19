#include "AmberFD.h"

AmberFD::AmberFD()
{
    n_sites = 0;
}

AmberFD::AmberFD(const int n_particles)
{
    nuclei.reserve(n_particles);
    frz_chg.reserve(n_particles);
    frz_exp.reserve(n_particles);
    dyn_exp.reserve(n_particles);
    pauli_coeff.reserve(n_particles);
    pauli_exp.reserve(n_particles);
    n_sites = 0;
}
AmberFD::~AmberFD()
{}

void AmberFD::add_particle(ParticleInfo parameters)
{
    nuclei.push_back(parameters.nuclei);
    frz_exp.push_back(parameters.frz_exp);
    frz_chg.push_back(parameters.frz_chg);
    dyn_exp.push_back(parameters.dyn_exp);
    pauli_exp.push_back(parameters.pauli_exp);
    pauli_coeff.push_back(parameters.pauli_coeff);
    n_sites += 1;
}

std::shared_ptr<FlucDens> AmberFD::create_fluc_dens_force()
{
    vec_d nuc(nuclei.begin(), nuclei.end());
    flucDens = std::shared_ptr<FlucDens>(new FlucDens(n_sites, &frz_chg[0], &nuc[0], &frz_exp[0], &dyn_exp[0]));
    return flucDens;
}

std::shared_ptr<DispersionPauli> AmberFD::create_disp_pauli_force()
{
    dispersionPauli = std::shared_ptr<DispersionPauli>(new DispersionPauli(n_sites, &nuclei[0], &pauli_exp[0], &pauli_coeff[0]));
    return dispersionPauli;
}