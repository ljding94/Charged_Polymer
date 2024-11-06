#ifndef _CHARGED_POLYMER_H
#define _CHARGED_POLYMER_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// hamiltonion parameters
struct Energy_parameter
{
    double kappa = 0; // bending modulii
    double A = 0;     // Yukawa strength: A
    double invK = 0;     // Yukawa range inverse: 1/K

    // Energy:
    // E = sum_i {\kappa/2*(t_i - t_{i-1})^2 + sum_{i,j} A*exp(-K(r_{ij}-1))/r_{ij}

};
struct observable
{
    // put all observable interested here
    double E;                            // energy + E_B + E_HSY
    double E_B;                           // total bending energy: bending energy
    double E_HSY;                         // total hard sphere + Yukawa energy
    double R2;                           // end to end distance square
    double Rg2;                           // radius of gyration

    std::vector<double> Sq{};  // structure factor
    std::vector<double> qB{};  // qB for Sq, B is the bead-bead distance, which is 1 in our unit
    std::vector<double> tts{}; // tangent-tangent correlation function, versus contour distance
    std::vector<double> spB{}; // s/B for tts calculation
};
struct bead
{
    // configuration related
    std::vector<double> r{0, 0, 0};     // position (x,y,z)
    std::vector<double> t{1.0, 0, 0.0}; // tangent, point to next bead R' = R+t, bead-bead distance B is unit vector
};

class charged_polymer
{
public:
    double beta;    // system temperature
    int L;          // length, (L+1) beads in total  (0)--(1)--(2)-- ... --(L-2)--(L-1)--(L)
    int fix_bead_0; // fix the 0th bead

    Energy_parameter Epar;
    std::vector<bead> polymer; // the polymer
    double E_sys;              // system energy
    double E_B_sys;             // system bending
    double E_HSY_sys;           // system Yukawa

    double d_theta; // for rot update

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_real_distribution<> rand_uni; // uniform distribution

    std::vector<double> rand_uni_vec(); // generate a random unit vector

    // initialization
    charged_polymer(double L_, Energy_parameter Epar_, bool random_Epar = 0);
    // MCMC update
    int update_bead_crankshaft(int bead_i, int bead_j);
    // con-rot update in a corn with angle d_theta
    int update_bead_pivot_right(int bead_i); // pivot right part
    // rotage the tangent of bead i, every bead > i position will be changed

    // check self-avoid condition
    int satisfy_self_avoiding_condition(int bead_i);                      // check if polymer[i] satisfy self-avoiding condition with all other beads
    int satisfy_self_avoiding_condition_by_group(int bead_i, int bead_j); //
    // check if polymer[:i] and polymer[i:] satisfy self-avoiding condition with all other beads
    double calc_Hard_Sphere_Yukawa_energy_by_group(int bead_i, int bead_j); //

    // some observable measurement
    observable measure_observable(int bin_num);
    // pair distribution function
    std::vector<double> calc_structure_factor(std::vector<double> qB);          // structure factor of the polymer, orientational averaged: ref: Pedersen 1996 equ 1
    std::vector<double> calc_tangent_pair_correlation(std::vector<double> spB); // calculate the pair tangent correlation distribution function
    double calc_radius_of_gyration_square();                                           // radius of gyration Rg^2

    // experiment
    void save_polymer_to_file(std::string filename);
    void save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble, bool save_detail = false);

    double run_MC_sweep(int step_per_sweep);
    void run_simultion(int therm_sweep, int MC_sweeps, int step_per_sweep, std::string folder, std::string finfo, int bin_num, int save_more_config);

    // spB stand for s per B or s/P

    // little tool
    std::vector<double> cross_product(std::vector<double> a, std::vector<double> b); // cross product of a and b
    double inner_product(std::vector<double> a, std::vector<double> b);              // inner product of a and b
    double calc_bending_of_two_t(std::vector<double> t1, std::vector<double> t2);    // calculate the bending angle between two tangent
    double calc_Hard_Sphere_Yukawa_of_two_r(std::vector<double> r1, std::vector<double> r2);     // calculate the Yukawa energy between two beads

    std::vector<double> Rodrigues_rotation(std::vector<double> v, std::vector<double> k, double theta); // v rotate around k by theta


};
#endif