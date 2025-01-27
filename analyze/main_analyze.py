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

    # versus invK, varying kappa
    if 1:
        folder = "../data/20250115"
        L = 500
        kappa = 5.0
        A = 5.0
        # invKs = np.concatenate((np.arange(1.0, 15.1, 1.0), np.arange(20.0, 100.1, 5.0)))
        invKs = np.arange(1.0, 10.1, 0.5)

        parameters = []
        finfos = []
        for kappa in [10, 20, 30, 40, 50, 60]:
            parameters = []
            finfos = []
            for invK in invKs:
                finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
                obs_filename = folder + f"/obs_{finfo}.csv"
                if not os.path.exists(obs_filename):
                    print(f"file {obs_filename} does not exist")
                    continue
                parameters.append((L, kappa, A, invK))
                finfos.append(finfo)
                # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
                # plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
                # plot_multi_config(folder, finfo)
                # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
            plot_obs(folder, finfos, parameters, "invK", f"_L{L}_kappa{kappa}_A{A}")
            # plot_obs(folder, finfos, parameters, "A", f"_L{L}_kappa{kappa}_invK{invK}")
        #return 0

    # versus invK, varying A
    if 0:
        folder = "../data/20250115"
        L = 500
        kappa = 30.0
        A = 5.0
        invKs = np.arange(1.5, 10.1, 0.5)
        parameters = []
        finfos = []
        for A in [1, 3, 5, 7, 9]:
            parameters = []
            finfos = []
            for invK in invKs:
                finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
                obs_filename = folder + f"/obs_{finfo}.csv"
                if not os.path.exists(obs_filename):
                    print(f"file {obs_filename} does not exist")
                    continue
                parameters.append((L, kappa, A, invK))
                finfos.append(finfo)
                # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
                # plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
                # plot_multi_config(folder, finfo)
                # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
            plot_obs(folder, finfos, parameters, "invK", f"_L{L}_kappa{kappa}_A{A}")
            # plot_obs(folder, finfos, parameters, "A", f"_L{L}_kappa{kappa}_invK{invK}")
        #return 0

    # versus invK, varying L
    if 1:
        folder = "../data/20250115"
        Ls = [300, 500, 700]
        kappa = 30.0
        A = 5.0
        # invKs = np.concatenate((np.arange(1.0, 15.1, 1.0), np.arange(20.0, 100.1, 5.0)))
        invKs = np.arange(1.5, 10.1, 0.5)

        parameters = []
        finfos = []
        for L in [300, 500, 700]:
            parameters = []
            finfos = []
            for invK in invKs:
                finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
                obs_filename = folder + f"/obs_{finfo}.csv"
                if not os.path.exists(obs_filename):
                    print(f"file {obs_filename} does not exist")
                    continue
                parameters.append((L, kappa, A, invK))
                finfos.append(finfo)
                # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
                # plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
                # plot_multi_config(folder, finfo)
                # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
            plot_obs(folder, finfos, parameters, "invK", f"_L{L}_kappa{kappa}_A{A}")
            # plot_obs(folder, finfos, parameters, "A", f"_L{L}_kappa{kappa}_invK{invK}")
        return 0




if __name__ == "__main__":
    main()
