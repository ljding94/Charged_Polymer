#!/opt/homebrew/bin/python3
from plot_analyze import *
from ML_analyze import *
import random
import time


def main():

    print("analyzing data using ML model")
    folder = "../data/20250116_rand"
    # folder = "../data/20241230_rand"
    rand_num = 20000
    # rand_max = 1000

    # filter by inK
    invK_ranges = [(i - 0.5, i + 0.5) for i in range(1, 10, 1)]
    # invK_ranges.append((3.25, 7.75))
    print("invK_ranges", invK_ranges)
    all_parameters_per_invK = [[] for _ in range(len(invK_ranges))]
    all_parameters = []
    gpr_A_invK_range = (3.5, 10.5)
    A_parameters = []
    gpr_kappa_invK_range = (0.5, 4.5)
    kappa_parameters = []
    for i in range(rand_num):
        filename = f"{folder}/obs_random_run{i}.csv"

        if os.path.exists(filename):
            data = np.genfromtxt(filename, delimiter=",", skip_header=1)
            L, kappa, A, invK = data[0, 1:5]

            for idx, (low, high) in enumerate(invK_ranges):
                if low <= invK < high:
                    all_parameters_per_invK[idx].append([i])
                    # break
            if gpr_A_invK_range[0] <= invK < gpr_A_invK_range[1]:
                A_parameters.append([i])
            if gpr_kappa_invK_range[0] <= invK < gpr_kappa_invK_range[1]:
                kappa_parameters.append([i])

            all_parameters.append([i])
        # if len(parameters) >= rand_max:
        #    break

    # print("parameters", parameters)
    print("total number of parameters", len(all_parameters))

    '''
    for idx, param_list in enumerate(all_parameters_per_invK):
        print(f"Total number of parameters in invK range {invK_ranges[idx]}: {len(param_list)}")

    # all data
    all_data_Vh = calc_svd(folder, all_parameters, note="_all", save_basis=True)

    for i in range(len(invK_ranges)):
        print("invK range", invK_ranges[i])
        calc_svd(folder, all_parameters_per_invK[i], use_Vh=all_data_Vh, note=f"_invK_{invK_ranges[i][0]}_{invK_ranges[i][1]}")
    invKs = [i for i in range(1, 10, 1)]
    NND_analysis(folder, invKs, all_parameters_per_invK, use_Vh=all_data_Vh)

    # plot_pddf_acf(folder, parameters, max_z=5, n_bin=100)
    #return 0
    '''

    # run gpr for two different invK ranges: all and gpr__A_invK_range

    # all param
    random.shuffle(all_parameters)
    all_parameters_train = all_parameters[: int(0.7 * len(all_parameters))]
    all_parameters_test = all_parameters[int(0.7 * len(all_parameters)) :]

    all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, all_parameters_train, ["R2", "Rg2", "kappa"], "_all")
    all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature = read_gp_and_feature_stats(folder, "_all")

    GaussianProcess_prediction(folder, all_parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature, ["R2", "Rg2", "kappa"], "_all")


    # A param
    random.shuffle(A_parameters)
    A_parameters_train = A_parameters[: int(0.7 * len(A_parameters))]
    A_parameters_test = A_parameters[int(0.7 * len(A_parameters)) :]
    print("np.shape(A_parameters)", np.shape(A_parameters))

    all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, A_parameters_train, ["A"], "_A")
    all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature = read_gp_and_feature_stats(folder, "_A")

    GaussianProcess_prediction(folder, A_parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature, ["A"], "_A")

    '''
    # kappa param
    random.shuffle(kappa_parameters)
    kappa_parameters_train = kappa_parameters[: int(0.7 * len(kappa_parameters))]
    kappa_parameters_test = kappa_parameters[int(0.7 * len(kappa_parameters)) :]
    print("np.shape(kappa_parameters)", np.shape(kappa_parameters))

    all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, kappa_parameters_train, ["kappa"], "_kappa")
    all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature = read_gp_and_feature_stats(folder, "_kappa")

    GaussianProcess_prediction(folder, kappa_parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature, ["kappa"], "_kappa")
    '''


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
