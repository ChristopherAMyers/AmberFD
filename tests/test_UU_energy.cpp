#include <stdio.h>
#include "FlucDens.h"
#include <vector>
#include "include/fileReader.h"
#include "include/asserts.h"
#include "include/common.h"
#include <map>
#include <utility>
#include <math.h>
#include <Vec3.h>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

std::map<std::string, int> symbol_to_nuc {{"H", 1}, {"C", 6}, {"N", 7}, {"O", 8}};
map_id nuc_to_exp_frz {{1, 2.50}, {6, 2.45}, {7, 2.17}, {8, 2.42}};
map_id nuc_to_exp_dyn {{1, 2.10}, {6, 2.05}, {7, 1.80}, {8, 1.90}};

void form_bonds(const vector<Vec3> coords, const vector<double> nuclei, vector<pair<int, int>> &bonds)
{
    bonds.clear();
    for (int i = 0; i < (int)coords.size(); i++)
    {
        for (int j = i+1; j < (int)coords.size(); j++)
        {
            bool add_bond = false;
            Vec3 dr = coords[j] - coords[i];
            double r = sqrt(dr.dot(dr));

            if ( ((int)nuclei[i] == 1) || ((int)nuclei[j] == 1) )
            {    if (r <= 1.2*ANG2BOHR)
                    add_bond = true;
            }
            else if (r <= 1.6*ANG2BOHR)
                add_bond = true;

            //printf(" %2d  %10.5f %2d  %2d  %d  %d\n", add_bond, r, i+1, j+1, (int)nuclei[i], (int)nuclei[j]);
            if (add_bond)
            {
                bonds.push_back(pair<int, int>(i, j));
                //printf(" bond %2d  %2d  %d  %d\n", i+1, j+1, (int)nuclei[i], (int)nuclei[j]);
            }
        }
    }
}

void print_exclusions(const FlucDens &fluc)
{
    for (int i = 0; i < fluc.get_num_frz_frz_exclusions(); i ++)
    {
        vector<int> ex_i;
        fluc.get_frz_frz_exclusions(i, ex_i);
        for (int j = 0; j < (int)ex_i.size(); j++)
            printf(" exclusion %3d  %3d \n", i + 1, ex_i[j] + 1);
    }
}


int main(int argc, char *argv[])
{
    double x, y, z;
    vec_d nuclei, frz_exp, dyn_exp, frz_chg, sites;
    vec_i frag1_idx, frag2_idx;
    vector<string> atoms;
    vector<Vec3> coords;
    int num_sites = 19;
    fileReader::flines_p data = fileReader::readFile("data/u_u_data_1.txt");
    int i = 0;
    for (auto line: data)
    {
        if (line[0][0] == '#') continue;
        nuclei.push_back(std::stod(line[1]));
        atoms.push_back(line[2]);
        x = std::stod(line[3]);
        y = std::stod(line[4]);
        z = std::stod(line[5]);
        coords.push_back(Vec3(x, y, z));
        sites.push_back(x);
        sites.push_back(y);
        sites.push_back(z);
        frz_chg.push_back(std::stod(line[6]));
        frz_exp.push_back(std::stod(line[7]));
        dyn_exp.push_back(std::stod(line[8]));
        if (i < num_sites) frag1_idx.push_back(i);
        else                frag2_idx.push_back(i);
        i++;
    }


    printf("There are %d sites\n", (int)nuclei.size());

    FlucDens fluc(atoms.size(), &frz_chg[0], &nuclei[0], &frz_exp[0], &dyn_exp[0]);
    fluc.set_dampening(1.5467, 1.4364);
    
    vector<pair<int, int>> bonds;
    form_bonds(coords, nuclei, bonds);
    fluc.create_frz_exclusions_from_bonds(bonds, 3);
    fluc.add_fragment(frag1_idx);
    fluc.add_fragment(frag2_idx);
    //print_exclusions(fluc);

    //fluc.print_params("frozen_exp parameters", "frozen_exp");
    //fluc.print_params("dynamic_exp parameters", "dynamic_exp");
    
    try{
        double energy = fluc.calc_energy(sites)*AU_2_KJ_PER_MOL;
        assert_equal_tol(energy, -163.9609916016171383, 1e-13);

        double pol_energy = fluc.calc_energy(sites, false)*AU_2_KJ_PER_MOL;
        assert_equal_tol(pol_energy, -10.0445370280187074, 1e-13);

        double frz_energy = fluc.calc_energy(sites, true, false)*AU_2_KJ_PER_MOL;
        assert_equal_tol(frz_energy, -153.9164545735984291, 1e-13);
        }
    catch(const std::exception& e) {
        cout << " exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
