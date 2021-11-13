#include <stdio.h>
#include "FlucDens.h"
#include <vector>
#include "include/fileReader.h"
#include "include/arrerts.h"
#include <map>
#include <utility>
#include <math.h>
#include <Vec3.h>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;
using std::pair;

std::map<std::string, int> symbol_to_nuc {{"H", 1}, {"C", 6}, {"N", 7}, {"O", 8}};
std::map<int, double> nuc_to_exp_frz {{1, 2.50}, {6, 2.45}, {7, 2.17}, {8, 2.42}};
std::map<int, double> nuc_to_exp_dyn {{1, 2.10}, {6, 2.05}, {7, 1.80}, {8, 1.90}};

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
            {    if (r <= 1.2)
                    add_bond = true;
            }
            else if (r <= 1.6)
                add_bond = true;

            //printf(" %2d  %10.5f %2d  %2d  %d  %d\n", add_bond, r, i+1, j+1, (int)nuclei[i], (int)nuclei[j]);
            if (add_bond)
            {
                bonds.push_back(pair<int, int>(i, j));
                printf(" bond %2d  %2d  %d  %d\n", i+1, j+1, (int)nuclei[i], (int)nuclei[j]);
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
    vector<Vec3> coords;
    vector<std::string> atoms;
    cout << "Loading molecule\n";
    fileReader::readXYZ("data/u_u.xyz", coords, atoms);
    cout << "Loading charges\n";
    fileReader::flines_p chg_data = fileReader::readFile("data/u_u.chg");
    
    
    cout << "Assigning data\n";
    vector<double> sites, nuclei, frz_exp, dyn_exp, frz_chg;
    for (int i = 0; i < (int)atoms.size(); i ++)
    {
        sites.push_back(coords[i][0]); sites.push_back(coords[i][1]); sites.push_back(coords[i][2]); 
        int nuc = symbol_to_nuc[atoms[i]];
        nuclei.push_back(nuc);
        frz_exp.push_back(nuc_to_exp_frz[nuc]);
        dyn_exp.push_back(nuc_to_exp_dyn[nuc]);
        frz_chg.push_back(std::stof(chg_data[i][0]));
    }
    // change exponents of the second molecule to create uneveness
    for (int i = (int)(atoms.size()/2); i < (int)atoms.size(); i++)
    {
        dyn_exp[i] += 0.25;
    }
    printf("There are %d atoms\n", (int)atoms.size());

    FlucDens fluc(atoms.size(), &frz_chg[0], &nuclei[0], &frz_exp[0], &dyn_exp[0]);
    fluc.set_dampening(1.5467, 1.4364);
    
    vector<pair<int, int>> bonds;
    form_bonds(coords, nuclei, bonds);
    fluc.create_frz_exclusions_from_bonds(bonds, 2);
    //print_exclusions(fluc);

    fluc.print_params("frozen_exp parameters", "frozen_exp");
    fluc.print_params("dynamic_exp parameters", "dynamic_exp");
    
    try{
        double energy = fluc.calc_energy(sites);
        assert_equal_tol(energy, -10.0619386256454657, 1e-13);

        double pol_energy = fluc.calc_energy(sites, false);
        assert_equal_tol(pol_energy, -5.517608011360994, 1e-13);

        double frz_energy = fluc.calc_energy(sites, true, false);
        assert_equal_tol(frz_energy, -4.5443306142844708, 1e-13);
        }
    catch(const std::exception& e) {
        cout << " exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
