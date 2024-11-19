#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np


def main():
    print("hello!")
    # plot_polymer_config("../data/20240828/config_L200_kappa10.0_f0.00_gL0.00.csv", "", True)
    # plot_R_distribution("../data/20240828/obs_MC_L200_kappa10.0_f0.00_gL0.00.csv")
    # plot_polymer_config("../data/scratch_local/20240829/config_L10_kappa10.0_f0.0_gL0.0.csv", "", True)
    # plot_R_distribution("../data/scratch_local/20240830/20240830/obs_MC_L200_kappa10.0_f0.00_gL0.00.csv")

    # plot_multi_config("../data/data_local/data_pool/L40_kappa5.0_A2.0_invK0.20")
    # return 0
    # folder = "../data/scratch_local/20240725"
    # plot_polymer_config(folder+"/config_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    # plot_polymer_config(folder+"/config_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # return 0

    folder = "../data/20241118"
    L = 500
    kappa = 5.0
    A = 1.0
    parameters = []
    finfos = []
    for kappa in [20, 50]:
        parameters = []
        finfos = []
        for invK in [1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0, 20.0, 40.0]:
            finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
            obs_filename = folder + f"/obs_{finfo}.csv"
            if not os.path.exists(obs_filename):
                continue
            parameters.append((L, kappa, A, invK))
            finfos.append(finfo)
            # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
            #plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
            # plot_multi_config(folder, finfo)
            # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
        plot_obs(folder, finfos, parameters, "invK", f"_L{L}_kappa{kappa}_A{A}")
    return 0

    folder = "../data/20240820"
    L = 200
    kappa = 10.0
    f = 0.0
    parameters = []
    finfos = []
    gLs = np.arange(0.00, 1.501, 0.10)
    for gL in gLs:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        parameters.append((L, kappa, f, gL))
        finfos.append(finfo)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
        # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
    # plot_obs(folder, finfos, parameters, "g")


if __name__ == "__main__":
    main()
