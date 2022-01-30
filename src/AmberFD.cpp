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

    omm_to_index[index] = n_sites - 1;
    index_to_omm[n_sites - 1] = index;

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
    DeltaR dR;
    if (periodicity.is_periodic)
        dR.getDeltaR(positions, (int)i*3, (int)j*3, periodicity);
    else
        dR.getDeltaR(positions, (int)i*3, (int)j*3);

    //  dispersion and pauli energies
    dispersionPauli->calc_one_pair(dR, i, j, pair_energies);

    //  fluctuating density and alectrostatics
    flucDens->calc_one_electro(dR, i, j, true, true, pair_energies);

    return pair_energies;
}

DeltaR AmberFD::getDeltaR(const vec_d &positions, int i, int j)
{
    if (periodicity.is_periodic)
        return DeltaR(positions, (int)i*3, (int)j*3, periodicity);
    else
        return DeltaR(positions, (int)i*3, (int)j*3);
}

void AmberFD::set_use_PBC(const bool is_periodic)
{
    periodicity.is_periodic = is_periodic;
    if (flucDens != nullptr)
        flucDens->set_use_PBC(is_periodic);
    if (dispersionPauli != nullptr)
        dispersionPauli->set_use_PBC(is_periodic);
}
void AmberFD::set_use_PBC(const bool is_periodic, const double x, const double y, const double z)
{
    //printf(" SETTING PBC 1: %d \n", is_periodic);
    periodicity.set(is_periodic, x, y, z);
    if (flucDens != nullptr)
        flucDens->set_use_PBC(is_periodic, x, y, z);
    if (dispersionPauli != nullptr)
        dispersionPauli->set_use_PBC(is_periodic, x, y, z);
    //printf(" SETTING PBC 2: %d \n", periodicity.is_periodic);
}

bool AmberFD::get_use_PBC()
{
    return periodicity.is_periodic;
}

Energies AmberFD::calc_energy_forces(const vec_d &positions)
{
    size_t i, j;
    total_energies.zero();
    Energies pair_energies;
    DeltaR dR;

    //  initialize solvers
    flucDens->initialize_calculation();
    dispersionPauli->initialize();
    
    for (i = 0; i < n_sites; i++)
    {
        for (j = i+1; j < n_sites; j++)
        {
            //  distances data
            if (periodicity.is_periodic)
                dR.getDeltaR(positions, (int)i*3, (int)j*3, periodicity);
            else
                dR.getDeltaR(positions, (int)i*3, (int)j*3);

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
        if (periodicity.is_periodic)
            flucDens->set_use_PBC(true, periodicity.box_size[0], periodicity.box_size[1], periodicity.box_size[2]);
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
        if (periodicity.is_periodic)
            dispersionPauli->set_use_PBC(true, periodicity.box_size[0], periodicity.box_size[1], periodicity.box_size[2]);
        return dispersionPauli;
    }
}

std::map<int, int> AmberFD::get_index_mapping()
{
    return omm_to_index;
}

void AmberFD::dump_to_file(std::string file_loc)
{
    FILE *file;
    file = fopen(file_loc.c_str(), "w");


    fprintf(file, "N_SITES %d\n", n_sites);
    double x = periodicity.box_size[0];
    double y = periodicity.box_size[1];
    double z = periodicity.box_size[2];
    fprintf(file, "PERIODICITY %d %.15e  %.15e  %.15e \n", periodicity.is_periodic, x, y, z);
    bool created_fluc = (flucDens != nullptr);
    bool created_disp = (dispersionPauli != nullptr);
    fprintf(file, "FLUC_FORCE %d\n", created_fluc);
    fprintf(file, "DISP_FORCE %d\n", created_disp);
    for(int i = 0; i < n_sites; i++)
    {
        int omm_index = index_to_omm[i];
        fprintf(file, "ATOM %d  %.15e  %.15e  %.15e  %.15e  %.15e  %d \n", nuclei[i], frz_chg[i], frz_exp[i], dyn_exp[i], pauli_exp[i], pauli_radii[i], omm_index);
    }
    fclose(file);
}

void AmberFD::load_from_file(std::string file_loc)
{
    using std::string;
    using std::stoi;
    using std::stod;
    string key, line;
    int num_sites = 0;
    bool created_fluc, created_disp;

    nuclei.clear();
    frz_chg.clear();
    frz_exp.clear();
    dyn_exp.clear();
    pauli_radii.clear();
    pauli_exp.clear();
    n_sites = 0;
    forces.clear();

    std::ifstream file;
    file.open(file_loc.c_str());
    int line_count = 0;
    std::vector<std::string> tokens;
    std::string token;
    if (file)
    {
        while(getline(file, line))
        {
                tokens.clear();
                std::istringstream ss(line);
                while(ss >> token)
                    tokens.push_back(token);
                key = tokens[0];
                if (key.compare("N_SITES") == 0)
                    num_sites == stoi(tokens[1]);
                else if (key.compare("PERIODICITY") == 0)
                {
                    set_use_PBC(stoi(tokens[1]), stod(tokens[2]), stod(tokens[3]), stod(tokens[4]));
                }
                else if (key.compare("ATOM") == 0)
                {
                    ParticleInfo particle(stoi(tokens[1]), stod(tokens[2]), stod(tokens[3]), stod(tokens[4]), stod(tokens[5]), stod(tokens[6]));

                    add_particle(stoi(tokens[7]), particle);
                    //printf("ADDING ATOM %d\n", n_sites);
                }
                else if (key.compare("FLUC_FORCE") == 0)
                    created_fluc = stoi(tokens[1]);
                else if (key.compare("DISP_FORCE") == 0)
                    created_disp = stoi(tokens[1]);
        }
    }

    if (created_fluc)
    {
        flucDens.reset();
        create_fluc_dens_force();
    }
    if (created_disp)
    {
        dispersionPauli.reset();
        create_disp_pauli_force();
    }

    file.close();
}