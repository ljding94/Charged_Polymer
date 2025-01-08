#!/opt/homebrew/bin/python3
from plot_analyze import *
from ML_analyze import *
import random
import time


def main():

    print("analyzing data using ML model")
    folder = "../data/20241230_rand"
    rand_num = 10000
    rand_max = 10000

    # filter by inK
    invK_ranges = [(i - 0.5, i + 0.5) for i in range(1, 17, 1)]
    #invK_ranges.append((3.25, 7.75))
    print("invK_ranges", invK_ranges)
    all_parameters_per_invK = [[] for _ in range(len(invK_ranges))]
    parameters = []
    for i in range(rand_num):
        filename = f"{folder}/obs_random_run{i}.csv"

        if os.path.exists(filename):
            data = np.genfromtxt(filename, delimiter=",", skip_header=1)
            L, kappa, A, invK = data[0, 1:5]

            for idx, (low, high) in enumerate(invK_ranges):
                if low <= invK < high:
                    all_parameters_per_invK[idx].append([i])
                    #break
            parameters.append([i])
        # if len(parameters) >= rand_max:
        #    break

    # print("parameters", parameters)
    print("total number of parameters", len(parameters))

    for idx, param_list in enumerate(all_parameters_per_invK):
        print(f"Total number of parameters in invK range {invK_ranges[idx]}: {len(param_list)}")

    # all data
    calc_svd(folder, parameters, note="_all")

    '''
    for i in range(len(invK_ranges)):
        print("invK range", invK_ranges[i])
        calc_svd(folder, all_parameters_per_invK[i], note=f"_invK_{invK_ranges[i][0]}_{invK_ranges[i][1]}")
    '''
    # plot_pddf_acf(folder, parameters, max_z=5, n_bin=100)
    return 0
    random.shuffle(parameters)
    parameters_train = parameters[: int(0.7 * len(parameters))]
    parameters_test = parameters[int(0.7 * len(parameters)) :]

    all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, parameters_train)
    all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature = read_gp_and_feature_stats(folder)

    GaussianProcess_prediction(folder, parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
