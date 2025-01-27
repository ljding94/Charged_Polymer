#include "charged_polymer.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <filesystem>
// #define PI 3.14159265358979323846
#define INF 1e9 // infinity

// initialization
charged_polymer::charged_polymer(double L_, Energy_parameter Epar_, bool random_Epar)
{
    // system related
    beta = 1;

    // set energy related
    // geometric

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_real_distribution<> rand_uni_set(0, 1);

    gen = gen_set;
    rand_uni = rand_uni_set;

    if (random_Epar)
    {
        L = 500 + int(0 * rand_uni(gen));
        Epar.kappa = 10 + 40 * rand_uni(gen);
        Epar.A = 10.0 * rand_uni(gen);
        Epar.invK = 1.0 + 9.0 * rand_uni(gen);
    }
    else
    {
        L = L_;
        Epar.kappa = Epar_.kappa;
        Epar.A = Epar_.A;
        Epar.invK = Epar_.invK;
    }

    d_theta = M_PI * 2.0 / 3 / (1.0 + std::sqrt(Epar.kappa));
    std::cout << "setting system param:" << "L" << L << ", kappa:" << Epar.kappa << ", A:" << Epar.A << ", 1/K:" << Epar.invK << "\n";

    // initial configuration
    polymer.resize(L + 1); // L+1 beads in total
    E_sys = 0;
    E_B_sys = 0;
    E_HSY_sys = 0;

    // initialize the polymer as straight along x axis,
    // such that there is no bending the shear energy
    // int mid_point = int(L / 2);
    for (int i = 0; i < L; i++)
    {
        polymer[i].r = {0, 0, 1.0 * i};
        polymer[i].t = {0, 0, 1.0};
    }
    polymer[L].r = {0, 0, 1.0 * L};
    polymer[L].t = {0, 0, 0}; // no tangent for the last bead

    // calc YKW energy
    for (int i = 0; i < L; i++)
    {
        for (int j = i + 1; j < L + 1; j++)
        {
            E_HSY_sys += calc_Hard_Sphere_Yukawa_of_two_r(polymer[i].r, polymer[j].r);
        }
    }
    if (E_HSY_sys > 0.5 * INF)
    {
        std::cout << "Error: E_HSY_sys > INF" << std::endl;
    }
    E_sys = E_B_sys + E_HSY_sys;

    std::cout << "Initial energies: E_sys = " << E_sys << ", E_B_sys = " << E_B_sys << ", E_HSY_sys = " << E_HSY_sys << std::endl;

    // above initial configuration give 0 bending and shear energy
}

int charged_polymer::update_bead_crankshaft(int bead_i, int bead_j)
{
    if (bead_i > bead_j)
    {
        // swap, make sure i is before j
        int bead = bead_i;
        bead_i = bead_j;
        bead_j = bead;
    }
    bead_j = std::min(bead_j, L);
    if (bead_j <= bead_i + 2)
    {
        return 0;
    }
    std::vector<bead> old_beads{}; // [bead_i, bead_j]
    double old_E = 0, new_E = 0;
    double old_E_B = 0, new_E_B = 0;
    double old_E_HSY = 0, new_E_HSY = 0;

    for (int j = bead_i; j <= bead_j; j++)
    {
        old_beads.push_back(polymer[j]); // (i),(i+1)... up to --(j),
    }
    // what updates:
    // t(i), r(i+1), t(i+1),..... r(j-1), t(j-1)
    // Energy change related:
    // bending (i-1)--(i)-- on (i) and  (j-1)--(j)-- on (i)
    // streching: no change since (0) and (L) are fixed
    // shearing (i) up to (j-1)

    old_E_B = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    old_E_B += (bead_j == L) ? 0 : calc_bending_of_two_t(polymer[bead_j - 1].t, polymer[bead_j].t);
    // bead_L has no t

    old_E_HSY = calc_Hard_Sphere_Yukawa_energy_by_group(bead_i, bead_j);
    if (old_E_HSY > 0.5 * INF)
    {
        std::cerr << "Warning: old_E_HSY > 0.5INF" << std::endl;
        return 0;
    }
    // TODO: consider storing all E_HSY_ij in a matrix, and update only the changed part to save time

    old_E = old_E_B + old_E_HSY;

    // update
    // 1. find the unit vector uni_t from (i) to (j)
    std::vector<double> uni_t = {polymer[bead_j].r[0] - polymer[bead_i].r[0],
                                 polymer[bead_j].r[1] - polymer[bead_i].r[1],
                                 polymer[bead_j].r[2] - polymer[bead_i].r[2]};

    double uni_t_norm = std::sqrt(uni_t[0] * uni_t[0] + uni_t[1] * uni_t[1] + uni_t[2] * uni_t[2]);
    uni_t[0] /= uni_t_norm;
    uni_t[1] /= uni_t_norm;
    uni_t[2] /= uni_t_norm;

    // 2. rotate everything between (i) and (j) about uni_t for random angle theta in renormalized [-d_theta,d_theta]
    double theta = d_theta / (2 * rand_uni(gen) - 1);
    // TODO: consider using simpler algorithm than Rodrigues_rotation

    std::vector<double> rik{0, 0, 0};
    double rik_norm;
    for (int k = bead_i; k < bead_j; k++)
    {
        // rik = {polymer[k].r[0] - polymer[bead_i].r[0], polymer[k].r[1] - polymer[bead_i].r[1], polymer[k].r[2] - polymer[bead_i].r[2]};
        // rik_norm = std::sqrt(rik[0] * rik[0] + rik[1] * rik[1] + rik[2] * rik[2]);
        polymer[k].t = Rodrigues_rotation(polymer[k].t, uni_t, theta); // this returned renomalized result
    }
    // reconstruct r
    for (int k = bead_i + 1; k < bead_j; k++)
    {
        polymer[k].r = {polymer[k - 1].r[0] + polymer[k - 1].t[0],
                        polymer[k - 1].r[1] + polymer[k - 1].t[1],
                        polymer[k - 1].r[2] + polymer[k - 1].t[2]};
    }
    if (std::abs(polymer[bead_j - 1].r[0] + polymer[bead_j - 1].t[0] - polymer[bead_j].r[0]) > 1e-6 ||
        std::abs(polymer[bead_j - 1].r[1] + polymer[bead_j - 1].t[1] - polymer[bead_j].r[1]) > 1e-6 ||
        std::abs(polymer[bead_j - 1].r[2] + polymer[bead_j - 1].t[2] - polymer[bead_j].r[2]) > 1e-6)
    {
        std::cout << "Error: r(j) not match after update bead_i" << std::endl;
        std::cout << "0:" << polymer[bead_j - 1].r[0] + polymer[bead_j - 1].t[0] << " " << polymer[bead_j].r[0] << std::endl;
        std::cout << "1:" << polymer[bead_j - 1].r[1] + polymer[bead_j - 1].t[1] << " " << polymer[bead_j].r[1] << std::endl;
        std::cout << "2:" << polymer[bead_j - 1].r[2] + polymer[bead_j - 1].t[2] << " " << polymer[bead_j].r[2] << std::endl;
        return 0;
    }
    // new energy
    new_E_B = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_E_B += (bead_j == L) ? 0 : calc_bending_of_two_t(polymer[bead_j - 1].t, polymer[bead_j].t);
    // bead_L has no t
    new_E_HSY = calc_Hard_Sphere_Yukawa_energy_by_group(bead_i, bead_j);
    new_E = new_E_B + new_E_HSY;

    // Metropolis
    if (new_E < 0.5 * INF && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        E_B_sys += new_E_B - old_E_B;
        E_HSY_sys += new_E_HSY - old_E_HSY;
        return 1;
    }
    else
    {
        for (int j = bead_i; j <= bead_j; j++)
        {
            polymer[j] = old_beads[j - bead_i];
        }

        return 0;
    }
}

int charged_polymer::update_bead_pivot_right(int bead_i)
{
    std::vector<bead> old_beads{}; // starting with bead_i
    double old_E = 0, new_E = 0;
    double old_E_B = 0, new_E_B = 0;
    double old_E_HSY = 0, new_E_HSY = 0;
    for (int j = bead_i; j < polymer.size(); j++)
    {
        old_beads.push_back(polymer[j]); // (i),(i+1)... up to (L-1),(L)
    }

    // for any bead_i with a t, which is [0,L-1]
    // everything j > bead_i us updated, since t_i rotate
    // what updates:
    // t(i), r(i+1), t(i+1), ..... r(L-1), t(L-1), r(L)
    // Energ change related
    // bending (i-1)--(i)--(i+1)-- on (i) if i>1 else, no bending changed
    // streching r(L) - r(i)
    // shearing all j>bead_i

    old_E_B = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    old_E_HSY = calc_Hard_Sphere_Yukawa_energy_by_group(bead_i, L + 1);
    if (old_E_HSY > 0.5 * INF)
    {
        std::cerr << "Warning: old_E_HSY > 0.5INF for pivot" << std::endl;
        return 0;
    }
    old_E = old_E_B + old_E_HSY;

    // update
    // 1. find random unit uni_v vector not parallel to t
    std::vector<double> uni_v = rand_uni_vec();
    while (std::abs(inner_product(uni_v, polymer[bead_i].t)) < 1e-6)
    { // uni_v must be not parallel to t
        uni_v = rand_uni_vec();
    }
    // 2. transform uni_v to vector uni_w perpendicular to t
    std::vector<double> uni_w = cross_product(polymer[bead_i].t, uni_v);
    double uni_w_norm = std::sqrt(uni_w[0] * uni_w[0] + uni_w[1] * uni_w[1] + uni_w[2] * uni_w[2]);
    uni_w[0] /= uni_w_norm;
    uni_w[1] /= uni_w_norm;
    uni_w[2] /= uni_w_norm;

    // also need uni_v perpendicular to both w and t
    uni_v = cross_product(polymer[bead_i].t, uni_w);
    double uni_v_norm = std::sqrt(uni_v[0] * uni_v[0] + uni_v[1] * uni_v[1] + uni_v[2] * uni_v[2]);
    uni_v[0] /= uni_v_norm;
    uni_v[1] /= uni_v_norm;
    uni_v[2] /= uni_v_norm;

    // 3. rotate every r(ij) towards uni_w for random angle theta such that cos(theta) uniform in [cos d_theta,1]
    // due to the length of rotated polymer, here use d_theta/(L-bead_i), to adjust acceptance rate
    double cos_theta = 1 - rand_uni(gen) * (1 - std::cos(1.0 * d_theta));
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    std::vector<double> r_ij(3, 0);
    double r_ij_norm = 0;
    double r_ij_norm_new = 0;
    double r_ij_v = 0;
    double r_ij_w = 0;
    double r_ij_t = 0;

    double tj_v = 0;
    double tj_w = 0;
    double tj_t = 0;
    double tj_norm = 0;
    // note: txw = v
    // check othorgonal of v,w,t
    if (std::abs(inner_product(uni_v, uni_w)) > 1e-6 || std::abs(inner_product(uni_v, polymer[bead_i].t)) > 1e-6 || std::abs(inner_product(uni_w, polymer[bead_i].t)) > 1e-6)
    {
        std::cout << "Error: uni_v, uni_w, t not orthogonal" << std::endl;
        std::cout << inner_product(uni_v, uni_w) << ":" << inner_product(uni_v, polymer[bead_i].t) << ":" << inner_product(uni_w, polymer[bead_i].t) << std::endl;
        return 0;
    }
    std::vector<double> ti_old = polymer[bead_i].t;
    for (int j = bead_i; j < polymer.size() - 1; j++) // polymer[L] has no t
    {
        // rotate tj towards uni_w by theta
        tj_v = inner_product(polymer[j].t, uni_v);
        tj_w = inner_product(polymer[j].t, uni_w);
        tj_t = inner_product(polymer[j].t, ti_old);
        polymer[j].t = {(cos_theta * tj_t - sin_theta * tj_w) * ti_old[0] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[0] + tj_v * uni_v[0],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[1] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[1] + tj_v * uni_v[1],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[2] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[2] + tj_v * uni_v[2]};
        // renormalize t
        tj_norm = std::sqrt(polymer[j].t[0] * polymer[j].t[0] + polymer[j].t[1] * polymer[j].t[1] + polymer[j].t[2] * polymer[j].t[2]);
        polymer[j].t[0] /= tj_norm;
        polymer[j].t[1] /= tj_norm;
        polymer[j].t[2] /= tj_norm;
    }

    // reconstruct r
    for (int k = bead_i + 1; k < polymer.size(); k++)
    {
        polymer[k].r = {polymer[k - 1].r[0] + polymer[k - 1].t[0],
                        polymer[k - 1].r[1] + polymer[k - 1].t[1],
                        polymer[k - 1].r[2] + polymer[k - 1].t[2]};
    }

    // calculate updated energy
    new_E_B = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_E_HSY = calc_Hard_Sphere_Yukawa_energy_by_group(bead_i, L + 1);
    new_E = new_E_B + new_E_HSY;

    if (new_E_HSY < 0.5 * INF && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        E_B_sys += new_E_B - old_E_B;
        E_HSY_sys += new_E_HSY - old_E_HSY;
        return 1;
    }
    else
    {
        for (int j = bead_i; j < polymer.size(); j++)
        {
            polymer[j] = old_beads[j - bead_i];
        }
        return 0;
    }
}

int charged_polymer::satisfy_self_avoiding_condition(int bead_i)
{
    double rij_norm2 = 0;
    std::vector<double> bead_r(3, 0);
    if (bead_i == L)
    {
        bead_r[0] = polymer[bead_i - 1].r[0] + polymer[bead_i - 1].t[0];
        bead_r[1] = polymer[bead_i - 1].r[1] + polymer[bead_i - 1].t[1];
        bead_r[2] = polymer[bead_i - 1].r[2] + polymer[bead_i - 1].t[2];
    }
    else
    {
        bead_r = polymer[bead_i].r;
    }
    // check distance between bead_i and all other beads
    for (int j = 0; j < L; j++)
    {
        if (j == bead_i - 1 || j == bead_i || j == bead_i + 1)
        {
            // neighbouring beads and self not need to check
            continue;
        }
        rij_norm2 = 0;
        for (int k = 0; k < 3; k++)
        {
            rij_norm2 += std::pow(bead_r[k] - polymer[j].r[k], 2);
        }
        // comparing the distance between i's next segment's center and j
        if (rij_norm2 < 1)
        {
            return 0;
        }
    }
    return 1;
}

int charged_polymer::satisfy_self_avoiding_condition_by_group(int bead_i, int bead_j)
{
    // check overlapping between j<bead_i and j>bead_i
    if (bead_i > bead_j)
    {
        // swap, make sure i is before j
        int bead = bead_i;
        bead_i = bead_j;
        bead_j = bead;
    }
    if (bead_j - bead_i < 2)
    {
        // nothing in between
        return 1;
    }
    double rjk_norm2 = 0;

    if (bead_j == L + 1)
    {
        for (int j = 0; j < bead_i; j++)
        {
            for (int k = bead_i + 1; k < polymer.size(); k++)
            {
                rjk_norm2 = 0;
                for (int l = 0; l < 3; l++)
                {
                    rjk_norm2 += std::pow(polymer[j].r[l] - polymer[k].r[l], 2);
                }
                if (rjk_norm2 < 1)
                {
                    return 0;
                }
            }
        }
    }
    else
    {
        for (int j = bead_i + 1; j < bead_j; j++)
        {
            for (int k = 0; k < polymer.size(); k++)
            {
                if (k > bead_i && k < bead_j)
                {
                    // within the group, no need to check
                    continue;
                }
                rjk_norm2 = 0;
                for (int l = 0; l < 3; l++)
                {
                    rjk_norm2 += (polymer[j].r[l] - polymer[k].r[l]) * (polymer[j].r[l] - polymer[k].r[l]);
                }
                if (rjk_norm2 < 1)
                {
                    return 0;
                }
            }
        }
    }
    return 1;
}

double charged_polymer::calc_Hard_Sphere_Yukawa_energy_by_group(int bead_i, int bead_j)
{
    // interaction between (bead_i,.., bead_j) and others: [0,...,bead_i] and, [bead_j,...,L]
    double E_HSY = 0;
    if (bead_j == L + 1)
    {
        for (int j = 0; j <= bead_i; j++)
        {
            for (int k = bead_i + 1; k < polymer.size(); k++)
            {
                // skip neighboring beads since their distance are fixed
                if (abs(j - k) == 1)
                {
                    continue;
                }
                E_HSY += calc_Hard_Sphere_Yukawa_of_two_r(polymer[j].r, polymer[k].r);
                if (E_HSY > 0.5 * INF)
                {
                    // std::cout << "j: " << j << ", k: " << k << ", polymer[j].r: {" << polymer[j].r[0] << ", " << polymer[j].r[1] << ", " << polymer[j].r[2] << "}, polymer[k].r: {" << polymer[k].r[0] << ", " << polymer[k].r[1] << ", " << polymer[k].r[2] << "}" << std::endl;
                    return INF;
                }
            }
        }
    }
    else
    {
        for (int j = bead_i + 1; j < bead_j; j++)
        {
            for (int k = 0; k < polymer.size(); k++)
            {
                if (k > bead_i && k < bead_j)
                {
                    // within the group, no need to check
                    continue;
                }
                if (abs(j - k) == 1)
                {
                    continue;
                }
                E_HSY += calc_Hard_Sphere_Yukawa_of_two_r(polymer[j].r, polymer[k].r);
                if (E_HSY > 0.5 * INF)
                {
                    return INF;
                }
            }
        }
    }
    return E_HSY;
}

void charged_polymer::save_polymer_to_file(std::string filename)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        f << "r[0],r[1],r[2],t[0],t[1],t[2]";

        for (int i = 0; i < size(polymer); i++)
        {
            f << "\n"
              << polymer[i].r[0] << "," << polymer[i].r[1] << "," << polymer[i].r[2]
              << "," << polymer[i].t[0] << "," << polymer[i].t[1] << "," << polymer[i].t[2];
        }
    }
    f.close();

} // end save_polymer_to_file

void charged_polymer::save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble, bool save_detail)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        int number_of_polymer = obs_ensemble.size();
        if (save_detail)
        {
            f << "beta=" << beta << "\n";
            f << "kappa=" << Epar.kappa << "\n";
            f << "A=" << Epar.A << "\n";
            f << "1/K=" << Epar.invK << "\n";
            f << "E(energy),E_B, E_HSY, R^2,Rg^2";
            for (int i = 0; i < number_of_polymer; i++)
            {
                f << "\n"
                  << obs_ensemble[i].E << "," << obs_ensemble[i].E_B << "," << obs_ensemble[i].E_HSY
                  << "," << obs_ensemble[i].R2 << "," << obs_ensemble[i].Rg2;
            }
        }
        else
        {
            // only save stats
            // find stats
            // calculate average and standard deviation of E
            double avg_E = 0.0, std_E = 0.0, avg_E_B = 0.0, std_E_B = 0.0, avg_E_HSY = 0.0, std_E_HSY = 0.0;
            double avg_R2 = 0.0, std_R2 = 0.0;
            double avg_Rg2 = 0.0, std_Rg2 = 0.0;

            std::vector<double> avg_Sq(obs_ensemble[0].Sq.size(), 0.0);
            std::vector<double> std_Sq(obs_ensemble[0].Sq.size(), 0.0);

            std::vector<double> avg_tts(obs_ensemble[0].tts.size(), 0.0);
            std::vector<double> std_tts(obs_ensemble[0].tts.size(), 0.0);

            double M = obs_ensemble.size();
            double sqrt_M = std::sqrt(M);

            // get the average
            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                avg_E += obs_ensemble[i].E;
                avg_E_B += obs_ensemble[i].E_B;
                avg_E_HSY += obs_ensemble[i].E_HSY;
                avg_R2 += obs_ensemble[i].R2;
                avg_R2 += obs_ensemble[i].R2;
                avg_Rg2 += obs_ensemble[i].Rg2;
                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    avg_Sq[j] += obs_ensemble[i].Sq[j];
                    avg_tts[j] += obs_ensemble[i].tts[j];
                }
            }
            avg_E /= M;
            avg_E_B /= M;
            avg_E_HSY /= M;
            avg_R2 /= M;
            avg_Rg2 /= M;

            if (obs_ensemble[0].Sq.size() != obs_ensemble[0].tts.size())
            {
                std::cout << "Error: Sq and tts size not match" << std::endl;
                return;
            }

            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                avg_Sq[j] /= M;
                avg_tts[j] /= M;
            }

            // get std
            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                std_E += (obs_ensemble[i].E - avg_E) * (obs_ensemble[i].E - avg_E);
                std_E_B += (obs_ensemble[i].E_B - avg_E_B) * (obs_ensemble[i].E_B - avg_E_B);
                std_E_HSY += (obs_ensemble[i].E_HSY - avg_E_HSY) * (obs_ensemble[i].E_HSY - avg_E_HSY);
                std_R2 += (obs_ensemble[i].R2 - avg_R2) * (obs_ensemble[i].R2 - avg_R2);
                std_Rg2 += (obs_ensemble[i].Rg2 - avg_Rg2) * (obs_ensemble[i].Rg2 - avg_Rg2);

                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    std_Sq[j] += (obs_ensemble[i].Sq[j] - avg_Sq[j]) * (obs_ensemble[i].Sq[j] - avg_Sq[j]);
                    std_tts[j] += (obs_ensemble[i].tts[j] - avg_tts[j]) * (obs_ensemble[i].tts[j] - avg_tts[j]);
                }
            }
            std_E = std::sqrt(std_E / M);
            std_E_B = std::sqrt(std_E_B / M);
            std_E_HSY = std::sqrt(std_E_HSY / M);
            std_R2 = std::sqrt(std_R2 / M);
            std_Rg2 = std::sqrt(std_Rg2 / M);

            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                std_Sq[j] = std::sqrt(std_Sq[j] / M);
                std_tts[j] = std::sqrt(std_tts[j] / M);
            }
            // write parameters and stats to the file
            f << "stats,L,kappa,A,1/K,E,E_B,E_HSY,R2,Rg2,Sq\n";
            f << "\nmean," << L << "," << Epar.kappa << "," << Epar.A << "," << Epar.invK << "," << avg_E
              << "," << avg_E_B << "," << avg_E_HSY << "," << avg_R2 << "," << avg_Rg2;
            for (int j = 0; j < avg_Sq.size(); j++)
            {
                f << "," << avg_Sq[j];
            }

            f << "\nstd/sqrt(number of polymer),NA,NA,NA,NA," << std_E / sqrt_M << "," << std_E_B / sqrt_M << "," << std_E_HSY / sqrt_M << "," << std_R2 / sqrt_M << "," << std_Rg2 / sqrt_M;
            for (int j = 0; j < std_Sq.size(); j++)
            {
                f << "," << std_Sq[j] / sqrt_M;
            }

            f << "\n qB,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].qB.size(); j++)
            {
                f << "," << obs_ensemble[0].qB[j];
            }

            f << "\n tts,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < avg_tts.size(); j++)
            {
                f << "," << avg_tts[j];
            }
            f << "\nstd_tts/sqrt(number of polymer),NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < std_tts.size(); j++)
            {
                f << "," << std_tts[j] / sqrt_M;
            }
            f << "\n r(for tts),NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].spB.size(); j++)
            {
                f << "," << obs_ensemble[0].spB[j];
            }
        }
    }
    f.close();
}

observable charged_polymer::measure_observable(int bin_num)
{
    // measure observable
    observable obs;
    obs.E = E_sys;
    obs.E_B = E_B_sys;
    obs.E_HSY = E_HSY_sys;

    std::vector<double> R = {0, 0, 0};
    R[0] = polymer[L - 1].r[0] + polymer[L - 1].t[0] - polymer[0].r[0];
    R[1] = polymer[L - 1].r[1] + polymer[L - 1].t[1] - polymer[0].r[1];
    R[2] = polymer[L - 1].r[2] + polymer[L - 1].t[2] - polymer[0].r[2];

    obs.R2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
    obs.Rg2 = calc_radius_of_gyration_square();

    std::vector<double> Sij = {0, 0, 0, 0, 0, 0};

    obs.qB.resize(bin_num);
    obs.Sq.resize(bin_num);

    double qB_i = 1e-3; // 0.1 / L;   // 0.2*M_PI/L; //0.1/L; ;
    double qB_f = 1;    // 100.0 / L; // M_PI;//100.0/L; //M_PI;
    for (int k = 0; k < bin_num; k++)
    {
        obs.qB[k] = qB_i * std::pow(qB_f / qB_i, 1.0 * k / (bin_num - 1)); // uniform in log scale
    }
    obs.Sq = calc_structure_factor(obs.qB);

    obs.spB.resize(bin_num);
    obs.tts.resize(bin_num);
    for (int k = 0; k < bin_num; k++)
    {
        obs.spB[k] = k; // 0, 1, 2, ..
    }
    obs.tts = calc_tangent_pair_correlation(obs.spB);

    return obs;
}

std::vector<double> charged_polymer::calc_structure_factor(std::vector<double> qB)
{
    // measure structure factor
    double q;
    int bin_num = qB.size();

    std::vector<std::vector<double>> R_all{}; // all scattering point's R[axis] value, including all scattering points in each segment

    for (int i = 0; i < polymer.size(); i++)
    {
        R_all.push_back(polymer[i].r);
    }
    int N = R_all.size();                      // total number of scattering points
    double r;                                  // for storating distance between two scattering points
    std::vector<double> SqB(bin_num, 1.0 / N); // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            r = std::sqrt(std::pow(R_all[i][0] - R_all[j][0], 2) + std::pow(R_all[i][1] - R_all[j][1], 2) + std::pow(R_all[i][2] - R_all[j][2], 2));
            // calculate S_q
            for (int k = 0; k < bin_num; k++)
            {
                q = qB[k];                                         // B = 1
                SqB[k] += 2.0 / N / N * std::sin(q * r) / (q * r); // normalize to 1 at q=0
            }
        }
    }
    return SqB;
}

std::vector<double> charged_polymer::calc_tangent_pair_correlation(std::vector<double> spB)
{
    int nbin = spB.size();
    std::vector<double> tts(nbin, 0);
    tts[0] = 1.0;                         // self correlation
    std::vector<int> pair_count(nbin, 0); // count number of t-t pair added to a s
    for (int i = 0; i < L - 1; i++)
    {
        for (int j = i + 1; j < L; j++)
        {
            if (j - i >= nbin)
            {
                continue;
            }
            tts[j - i] += inner_product(polymer[i].t, polymer[j].t);
            pair_count[j - i]++;
        }
    }
    for (int i = 1; i < nbin; i++)
    {
        tts[i] /= pair_count[i];
    }
    return tts;
}

double charged_polymer::calc_radius_of_gyration_square()
{
    // calculate radius of gyration
    double Rg2 = 0;
    // find center of mass
    double x_c = 0, y_c = 0, z_c = 0;
    for (int i = 0; i < polymer.size(); i++)
    {
        x_c = x_c + polymer[i].r[0];
        y_c = y_c + polymer[i].r[1];
        z_c = z_c + polymer[i].r[2];
    }
    x_c = x_c / L;
    y_c = y_c / L;
    z_c = z_c / L;

    // calculate Rg
    for (int i = 0; i < polymer.size(); i++)
    {
        Rg2 += std::pow(polymer[i].r[0] - x_c, 2) + std::pow(polymer[i].r[1] - y_c, 2) + std::pow(polymer[i].r[2] - z_c, 2);
    }
    Rg2 /= polymer.size();
    return Rg2;
}

std::vector<double> charged_polymer::cross_product(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> c(3, 0);
    for (int j = 0; j < 3; j++)
    {
        c[j] = a[(j + 1) % 3] * b[(j + 2) % 3] - a[(j + 2) % 3] * b[(j + 1) % 3];
    }
    return c;
}

double charged_polymer::inner_product(std::vector<double> a, std::vector<double> b)
{
    double c = 0;
    for (int j = 0; j < 3; j++)
    {
        c += a[j] * b[j];
    }
    return c;
}

std::vector<double> charged_polymer::Rodrigues_rotation(std::vector<double> v, std::vector<double> k, double theta)
{
    // v rotate around k by theta
    std::vector<double> v_rot(3, 0);
    // 1. find (k x v), k cross v and kv, k.v
    std::vector<double> kxv = cross_product(k, v);
    double kv = inner_product(k, v);
    // 2. do the calculation
    for (int j = 0; j < 3; j++)
    {
        v_rot[j] = v[j] * std::cos(theta) + kxv[j] * std::sin(theta) + k[j] * kv * (1 - std::cos(theta));
    }
    // 3. renormalize v_rot
    double v_rot_norm = std::sqrt(v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2]);
    v_rot[0] /= v_rot_norm;
    v_rot[1] /= v_rot_norm;
    v_rot[2] /= v_rot_norm;

    if (std::abs(v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2] - 1) > 1e-6)
    {
        printf("v_rot^2=%f\n", v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2]);
        throw std::runtime_error("Error: v_rot^2 is not 1.");
    }

    return v_rot;
}

std::vector<double> charged_polymer::rand_uni_vec()
{
    // 1. generat random point in unit ball
    double x, y, z;
    do
    {
        x = 2 * rand_uni(gen) - 1;
        y = 2 * rand_uni(gen) - 1;
        z = 2 * rand_uni(gen) - 1;
    } while (x * x + y * y + z * z > 1);
    // 2. normalize
    double norm = std::sqrt(x * x + y * y + z * z);
    x /= norm;
    y /= norm;
    z /= norm;
    return {x, y, z};
}

double charged_polymer::calc_bending_of_two_t(std::vector<double> t1, std::vector<double> t2)
{
    double bend = 0;
    bend += (t1[0] - t2[0]) * (t1[0] - t2[0]);
    bend += (t1[1] - t2[1]) * (t1[1] - t2[1]);
    bend += (t1[2] - t2[2]) * (t1[2] - t2[2]);
    return 0.5 * Epar.kappa * bend;
}
double charged_polymer::calc_Hard_Sphere_Yukawa_of_two_r(std::vector<double> r1, std::vector<double> r2)
{
    double r = 0;
    r = (r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]);
    r = std::sqrt(r);
    if (r < 1)
    {
        return INF;
    }
    if (Epar.invK == 0 || r > 5.0 * Epar.invK)
    {
        return 0;
    }

    return Epar.A * std::exp(-r / Epar.invK) / r;
}

double charged_polymer::run_MC_sweep(int step_per_sweep)
{
    int bead_i, bead_j, bead_ij;
    double acceptance_rate = 0;
    // int mid_point = int(L / 2);
    for (int n = 0; n < step_per_sweep; n++)
    {
        bead_ij = 2 + int(rand_uni(gen) * Epar.kappa); //~[2,L]
        bead_i = int(rand_uni(gen) * (L - bead_ij));
        bead_i += bead_ij;
        bead_j = bead_i + bead_ij;

        acceptance_rate += update_bead_crankshaft(bead_i, bead_j);
        /*
        // take between [0,mid_poin-bead_ij] and [mid_point,L-bead_ij]
        if (bead_i > mid_point - bead_ij)
        {
            bead_i += bead_ij;
            bead_j = bead_i + bead_ij;
            bead_j = std::min(bead_j, L - 1);
        }
        else
        {
            bead_j = bead_i + bead_ij;
            bead_j = std::min(bead_j, mid_point);
        }
        */

        bead_i = int(rand_uni(gen) * L); // take from [1,L-1]
        update_bead_pivot_right(bead_i);
    }
    return acceptance_rate / step_per_sweep;
}

void charged_polymer::run_simultion(int therm_sweep, int MC_sweeps, int step_per_sweep, std::string folder, std::string finfo, int bin_num, int save_more_config)
{
    std::vector<observable> obs_ensemble;

    double acceptance_rate = 0;

    beta = 0; // randomization
    for (int i = 0; i < therm_sweep; i++)
    {
        std::cout << "randomization(beta=0) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        run_MC_sweep(step_per_sweep);
    }

    std::cout << "\n";
    // thermalization using tempering
    for (int i = 0; i < therm_sweep; i++)
    {
        beta = (i + 1) * 1.0 / therm_sweep; // tempering
        std::cout << "thermalization/tempering (beta -> 1) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        run_MC_sweep(step_per_sweep);
    }

    // data collection run
    std::cout << "\n";
    if (beta != 1)
    {
        std::cout << "Error: beta is not 1 at the end of thermalization\n";
        beta = 1;
    }
    if (save_more_config)
    {
        std::cout << "saving more config, creating foler: " << folder + "/" + finfo << std::endl;
        std::filesystem::create_directory(folder + "/" + finfo);
    }

    for (int i = 0; i < MC_sweeps; i++)
    {
        std::cout << "MC sweep " << i << " out of " << MC_sweeps << " (" << (i * 100) / MC_sweeps << "%)\r";
        acceptance_rate += run_MC_sweep(step_per_sweep);

        obs_ensemble.push_back(measure_observable(bin_num));
        if (i % 100 == 0 && save_more_config)
        {
            save_polymer_to_file(folder + "/" + finfo + "/config_" + std::to_string(int(i / 100)) + ".csv");
        }
    }
    std::cout << "\n";
    std::cout << "crankshaft_acceptance rate:" << acceptance_rate / MC_sweeps << std::endl;

    save_polymer_to_file(folder + "/config_" + finfo + ".csv");
    // save_observable_to_file(folder + "/obs_MC_" + finfo + ".csv", obs_ensemble, true);
    save_observable_to_file(folder + "/obs_" + finfo + ".csv", obs_ensemble, false);
    // save_observable_to_file(folder + "/obs_MC_" + finfo + ".csv", obs_ensemble, true);
}