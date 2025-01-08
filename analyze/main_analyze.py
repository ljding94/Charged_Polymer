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

    # versus invK
    if 0:
        folder = "../data/20241214"
        L = 300
        kappa = 5.0
        A = 5.0
        # invKs = np.concatenate((np.arange(1.0, 15.1, 1.0), np.arange(20.0, 100.1, 5.0)))
        invKs = np.arange(1.0, 20.1, 1.0)

        parameters = []
        finfos = []
        for kappa in [5, 10, 15, 20, 25, 30]:
            parameters = []
            finfos = []
            for invK in invKs:
                finfo = f"L{L}_kappa{kappa:.0f}_A{A:.1f}_invK{invK:.1f}"
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

    # versus A
    if 1:
        folder = "../data/20241214"
        L = 300
        kappa = 5.0
        As = np.arange(1.0, 20.1, 1.0)
        invK = 10.0
        parameters = []
        finfos = []
        for kappa in [5, 10, 15, 20, 25, 30]:
            for invK in [2, 4, 6, 8]:
                parameters = []
                finfos = []
                for A in As:
                    finfo = f"L{L}_kappa{kappa:.0f}_A{A:.1f}_invK{invK:.1f}"
                    obs_filename = folder + f"/obs_{finfo}.csv"
                    if not os.path.exists(obs_filename):
                        continue
                    parameters.append((L, kappa, A, invK))
                    finfos.append(finfo)
                    # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
                    # plot_polymer_config(folder + f"/config_{finfo}.csv", finfo)
                    # plot_multi_config(folder, finfo)
                    # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
                plot_obs(folder, finfos, parameters, "A", f"_L{L}_kappa{kappa}_invK{invK}")
                # plot_obs(folder, finfos, parameters, "A", f"_L{L}_kappa{kappa}_invK{invK}")
        return 0


if __name__ == "__main__":
    main()
