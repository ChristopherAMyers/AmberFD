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

    forces.reserve(n_particles);
}
AmberFD::~AmberFD()
{}

int AmberFD::add_particle(int index, ParticleInfo parameters)
{
    nuclei.push_back(parameters.nuclei);
    frz_exp.push_back(parameters.frz_exp);
    frz_chg.push_back(parameters.frz_chg);
    dyn_exp.push_back(parameters.dyn_exp);
    pauli_exp.push_back(parameters.pauli_exp);
    pauli_radii.push_back(parameters.pauli_radii);
    forces.push_back(vec_d(3, 0.0));
    n_sites += 1;

    system_idx[index] = n_sites - 1;

    return n_sites - 1;
}

void AmberFD::add_fragment(const vec_i frag_idx)
{
    dispersionPauli->create_exclusions_from_fragment(frag_idx);
    flucDens->add_fragment(frag_idx);
}

std::vector<vec_d> AmberFD::get_forces()
{
    return forces;
}

Energies AmberFD::calc_one_pair(const vec_d &positions, int i, int j)
{
    Energies pair_energies;
    DeltaR dR(positions, (int)i*3, (int)j*3);

    //  dispersion and pauli energies
    dispersionPauli->calc_one_pair(dR, i, j, pair_energies);

    //  fluctuating density and alectrostatics
    flucDens->calc_one_electro(dR, i, j, true, true, pair_energies);

    return pair_energies;
}

Energies AmberFD::calc_energy_forces(const vec_d &positions)
{
    size_t i, j;
    total_energies.zero();
    Energies pair_energies;

    //  initialize solvers
    flucDens->initialize_calculation();
    dispersionPauli->initialize();
    
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances data
            DeltaR dR(positions, (int)i*3, (int)j*3);

            //  dispersion and pauli energies
            dispersionPauli->calc_one_pair(dR, i, j, pair_energies);
            total_energies.pauli += pair_energies.pauli;
            total_energies.disp += pair_energies.disp;
            total_energies.pauli_wall += pair_energies.pauli_wall;

            //  fluctuating density and alectrostatics
            flucDens->calc_one_electro(dR, i, j, true, true, pair_energies);
            total_energies.elec_elec += pair_energies.elec_elec;
            total_energies.elec_nuc += pair_energies.elec_nuc;
            total_energies.nuc_nuc += pair_energies.nuc_nuc;
            total_energies.frz += pair_energies.frz;
            total_energies.vct += pair_energies.vct;
        }
    }
    //  minimize fluc-dens energy
    flucDens->solve_minimization();
    total_energies.pol = flucDens->get_polarization_energy();
    total_energies.vct = flucDens->get_ct_energy();

    //  copy over forces from the solvers
    std::vector<vec_d> fluc_forces = flucDens->get_forces();
    std::vector<vec_d> disp_forces = dispersionPauli->get_forces();
    for(size_t i = 0; i < forces.size(); i++)
    {
        for(int j = 0; j < 3; j++)
            forces[i][j] = fluc_forces[i][j] + disp_forces[i][j];
    }

    return total_energies;
}

std::shared_ptr<FlucDens> AmberFD::get_fluc_dens_force(bool create_if_null)
{
    if (flucDens == nullptr)
        if (create_if_null)
            return create_fluc_dens_force();
        else
            throw std::runtime_error("FlucForce has not been created");
    else
        return flucDens;
}

std::shared_ptr<DispersionPauli> AmberFD::get_disp_pauli_force(bool create_if_null)
{
    if (dispersionPauli == nullptr)
        if (create_if_null)
            return create_disp_pauli_force();
        else
            throw std::runtime_error("DispersionPauli has not been created");
    else
        return dispersionPauli;
}

std::shared_ptr<FlucDens> AmberFD::create_fluc_dens_force()
{
    //  create a new
    if (flucDens != nullptr)
    {
        throw std::runtime_error("FlucForce has already been created");
    }
    else
    {
        vec_d nuc(nuclei.begin(), nuclei.end());
        flucDens = std::shared_ptr<FlucDens>(new FlucDens(n_sites, &frz_chg[0], &nuc[0], &frz_exp[0], &dyn_exp[0]));
        return flucDens;
    }
}

std::shared_ptr<DispersionPauli> AmberFD::create_disp_pauli_force()
{
    if (dispersionPauli != nullptr)
    {
        throw std::runtime_error("DispersionPauli has already been created");
    }
    else
    {
        dispersionPauli = std::shared_ptr<DispersionPauli>(new DispersionPauli(n_sites, &nuclei[0], &pauli_exp[0], &pauli_radii[0]));
        return dispersionPauli;
    }
}

std::map<int, int> AmberFD::get_index_mapping()
{
    return system_idx;
}