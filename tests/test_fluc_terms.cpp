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

            if (add_bond)
            {
                bonds.push_back(pair<int, int>(i, j));
            }
        }
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

    FlucDens fluc(atoms.size(), &frz_chg[0], &nuclei[0], &frz_exp[0], &dyn_exp[0]);
    fluc.set_dampening(1.5467, 1.4364);
    
    vector<pair<int, int>> bonds;
    form_bonds(coords, nuclei, bonds);
    fluc.create_frz_exclusions_from_bonds(bonds, 3);
    fluc.add_fragment(frag1_idx);
    fluc.add_fragment(frag2_idx);

    double inv_r, a, b, exp_ar, exp_br, dEdR;
    inv_r = 1.0/3.0;
    a = 2.0;
    b = 2.5;
    exp_ar = exp(-a/inv_r);
    exp_br = exp(-b/inv_r);
    
    double ee_1 = fluc.elec_elec_energy(inv_r, a, b, exp_ar, exp_br, dEdR);
    double eZ_1a = fluc.elec_nuclei_energy(inv_r, a, exp_ar, dEdR);
    double eZ_1b = fluc.elec_nuclei_energy(inv_r, b, exp_br, dEdR);
    
    inv_r = 1/5.0;
    exp_ar = exp(-a/inv_r);
    exp_br = exp(-b/inv_r);
    double ee_2 = fluc.elec_elec_energy(inv_r, a, b, exp_ar, exp_br, dEdR);
    double eZ_2a = fluc.elec_nuclei_energy(inv_r, a, exp_ar, dEdR);
    double eZ_2b = fluc.elec_nuclei_energy(inv_r, b, exp_br, dEdR);

    b = 2.0;
    exp_br = exp(-b/inv_r);
    double ee_3 = fluc.elec_elec_energy(inv_r, a, b, exp_ar, exp_br, dEdR);
    double eZ_3a = fluc.elec_nuclei_energy(inv_r, a, exp_ar, dEdR);
    double eZ_3b = fluc.elec_nuclei_energy(inv_r, b, exp_br, dEdR);

    a = 1.0;
    exp_ar = exp(-a/inv_r);
    double ee_4 = fluc.elec_elec_energy(inv_r, a, b, exp_ar, exp_br, dEdR);
    double eZ_4a = fluc.elec_nuclei_energy(inv_r, a, exp_ar, dEdR);
    double eZ_4b = fluc.elec_nuclei_energy(inv_r, b, exp_br, dEdR);

    a = 1.0;
    b = a + 0.0009;
    exp_ar = exp(-a/inv_r);
    exp_br = exp(-b/inv_r);
    double ee_5 = fluc.elec_elec_energy(inv_r, a, b, exp_ar, exp_br, dEdR);
    double eZ_5a = fluc.elec_nuclei_energy(inv_r, a, exp_ar, dEdR);
    double eZ_5b = fluc.elec_nuclei_energy(inv_r, b, exp_br, dEdR);

    
    printf(" ee_1  %20.16f \n", ee_1);
    printf(" eZ_1a %20.16f \n", eZ_1a);
    printf(" eZ_1b %20.16f \n", eZ_1b);

    printf(" ee_2  %20.16f \n", ee_2);
    printf(" eZ_2a %20.16f \n", eZ_2a);
    printf(" eZ_2b %20.16f \n", eZ_2b);

    printf(" ee_3  %20.16f \n", ee_3);
    printf(" eZ_3a %20.16f \n", eZ_3a);
    printf(" eZ_3b %20.16f \n", eZ_3b);

    printf(" ee_4  %20.16f \n", ee_4);
    printf(" eZ_4a %20.16f \n", eZ_4a);
    printf(" eZ_4b %20.16f \n", eZ_4b);

    printf(" ee_5  %20.16f \n", ee_5);
    printf(" eZ_5a %20.16f \n", eZ_5a);
    printf(" eZ_5b %20.16f \n", eZ_5b);

    double tol = 1e-13;
    try{
        assert_equal_tol(ee_1    ,   0.3244950545375144, tol);
        assert_equal_tol(eZ_1a   ,  -0.3300283304311115, tol);
        assert_equal_tol(eZ_1b   ,  -0.3324576164139326, tol);
        assert_equal_tol(ee_2    ,   0.1997985733307481, tol);
        assert_equal_tol(eZ_2a   ,  -0.1999455200842850, tol);
        assert_equal_tol(eZ_2b   ,  -0.1999945963529005, tol);
        assert_equal_tol(ee_3    ,   0.1995690790000044, tol);
        assert_equal_tol(eZ_3a   ,  -0.1999455200842850, tol);
        assert_equal_tol(eZ_3b   ,  -0.1999455200842850, tol);
        assert_equal_tol(ee_4    ,   0.1932033985922263, tol);
        assert_equal_tol(eZ_4a   ,  -0.1952834371006402, tol);
        assert_equal_tol(eZ_4b   ,  -0.1999455200842850, tol);
        assert_equal_tol(ee_5    ,   0.1842169249586514, tol);
        assert_equal_tol(eZ_5a   ,  -0.1952834371006402, tol);
        assert_equal_tol(eZ_5b   ,  -0.1953015954875795, tol);
    }
    catch(const std::exception e){
        cout << " exception: " << e.what() << endl;
        return 1;
    }

    return 0;


}
