// Copyright[2024] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>
#include "charged_polymer.h"

/*
plan:
for various L:
1. effect of A,K, on persistence length
2. effect of A,K, on Rg2
3. effect of A,K, on R2
4. effect of A,K, on Sq
*/

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();

    double beta = 1;
    Energy_parameter Epar;

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm *timeinfo = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    std::string today(buffer);
    std::cout << today << std::endl;

    // precision run with specified parameters
    int L;
    int n_index = 0;
    int save_more_config = 0;
    std::string folder = "./";
    std::string finfo;
    if (argc == 7)
    {
        std::cout << "Usage: " << argv[0] << " L kappa f g save_more_config folder" << std::endl;

        L = std::atoi(argv[1]);
        Epar.kappa = std::atof(argv[2]);
        Epar.A = std::atof(argv[3]);
        Epar.invK = std::atof(argv[4]);
        save_more_config = std::atoi(argv[5]);
        folder = argv[6];

        charged_polymer polymer(L, Epar, 0);
        finfo = "L" + std::string(argv[1]) + "_kappa" + std::string(argv[2]) + "_A" + std::string(argv[3]) + "_invK" + std::string(argv[4]);
        int bin_num;
        int therm_sweeps;
        int MC_sweeps;
        int step_per_sweep;

        bin_num = 100;
        therm_sweeps = 2000;      // 1500
        MC_sweeps = 4000;        // 3000
        step_per_sweep = L; // 100 * L;

        // polymer.save_polymer_to_file(folder + "/config_" + finfo + "_init.csv");
        polymer.save_polymer_to_file(folder + "/config_init_" + finfo + ".csv"); // save sample polymer
        polymer.run_simultion(therm_sweeps, MC_sweeps, step_per_sweep, folder, finfo, bin_num, save_more_config);
    }
    else if (argc == 5)
    {
        L = std::atoi(argv[1]);
        n_index = std::atoi(argv[2]);
        save_more_config = std::atoi(argv[3]);
        folder = argv[4];
        charged_polymer polymer(L, Epar, 1);
        finfo = "L" + std::string(argv[1]) + "_random_run" + std::string(argv[2]);

        int bin_num;
        int therm_sweeps;
        int MC_sweeps;
        int step_per_sweep;

        bin_num = 100;
        therm_sweeps = 2000;
        MC_sweeps = 4000;
        step_per_sweep = L;

        // polymer.save_polymer_to_file(folder + "/config_" + finfo + "_init.csv");
        polymer.save_polymer_to_file(folder + "/config_init_" + finfo + ".csv"); // save sample polymer
        polymer.run_simultion(therm_sweeps, MC_sweeps, step_per_sweep, folder, finfo, bin_num, save_more_config);
    }
    else
    {
        std::cout << "input error\n";
        return 0;
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << finfo << ": Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}