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
    pauli_radii.reserve(n_particles);
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
    pauli_radii.push_back(parameters.pauli_radii);
    n_sites += 1;
}

void AmberFD::add_fragment(const vec_i frag_idx)
{
    dispersionPauli->create_exclusions_from_fragment(frag_idx);
    flucDens->add_fragment(frag_idx);
}

Energies AmberFD::calc_energy_forces(const vec_d &positions)
{
    size_t i, j;
    total_energies.reset();
    Energies pair_energies;

    //  initialize solvers
    flucDens->initialize_calculation();
    dispersionPauli->initialize();
    
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances data
            DeltaR dR(positions, (int)j*3, (int)i*3);

            //  dispersion and pauli energies
            dispersionPauli->calc_one_pair(dR, i, j, pair_energies);
            total_energies.pauli += pair_energies.pauli;
            total_energies.disp += pair_energies.disp;

            //  fluctuating density and alectrostatics
            flucDens->calc_one_electro(dR, i, j, true, true, pair_energies);
            total_energies.frz += pair_energies.frz;
            total_energies.vct += pair_energies.vct;
        }
    }
    //  minimize fluc-dens energy
    flucDens->solve_minimization();
    total_energies.pol = flucDens->get_polarization_energy();

    return total_energies;
}

std::shared_ptr<FlucDens> AmberFD::create_fluc_dens_force()
{
    vec_d nuc(nuclei.begin(), nuclei.end());
    flucDens = std::shared_ptr<FlucDens>(new FlucDens(n_sites, &frz_chg[0], &nuc[0], &frz_exp[0], &dyn_exp[0]));
    return flucDens;
}

std::shared_ptr<DispersionPauli> AmberFD::create_disp_pauli_force()
{
    dispersionPauli = std::shared_ptr<DispersionPauli>(new DispersionPauli(n_sites, &nuclei[0], &pauli_exp[0], &pauli_radii[0]));
    return dispersionPauli;
}