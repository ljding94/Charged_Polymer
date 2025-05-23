import numpy as np
from plot_analyze import *
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import os
from scipy.optimize import curve_fit
import pickle
from scipy.spatial import cKDTree
from plot_analyze import calc_Sq_discrete_infinite_thin_rod


def get_all_feature_Sq_data(folder, parameters):

    all_L, all_kappa, all_A, all_invK = [], [], [], []
    all_R2, all_Rg2 = [], []  # R2, Rg2 related
    all_Sq = []
    qB = []
    for i in range(len(parameters)):
        run_num = parameters[i][0]
        filename = f"{folder}/obs_random_run{run_num}.csv"
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        L, kappa, A, invK = data[0, 1:5]
        R2, Rg2 = data[0, 8:10]

        all_L.append(L)
        all_kappa.append(kappa)
        all_A.append(A)
        all_invK.append(invK)
        all_R2.append(R2 / L**2)
        all_Rg2.append(Rg2 / L**2)

        Sq = data[0, 10:]
        all_Sq.append(Sq)
        qB = data[2, 10:]

    Sq_rod = calc_Sq_discrete_infinite_thin_rod(qB, int(L))
    all_feature = np.array([all_L, all_kappa, all_A, all_invK, all_R2, all_Rg2])
    all_feature_name = ["L", "kappa", "A", "invK", "R2", "Rg2"]
    all_feature_tex = [r"L", r"$\kappa$", r"$A$", r"$1/K$", r"$R^2/L^2$", r"$R_g^2/L^2$"]
    all_Sq = np.array(all_Sq)
    all_Sq = np.log(all_Sq)
    #all_Sq = all_Sq/np.array(Sq_rod)
    qB = np.array(qB)
    return all_feature.T, all_feature_name, all_feature_tex, all_Sq, qB


def calc_nearest_neighbor_distance(SqV, C):

    # Step 1: Build a k-d tree for efficient neighbor search
    tree = cKDTree(SqV)  # Use only spatial coordinates (x, y, z)
    distances, indices = tree.query(SqV, k=2)  # Find nearest neighbors (k=2)

    # Step 2: Compute color differences
    color_differences = np.abs(C - C[indices[:, 1]])

    # Step 3: Normalize by color range
    color_min = np.min(C)  # Minimum color value
    color_max = np.max(C)  # Maximum color value
    color_range = color_max - color_min  # Range of color values

    # Avoid division by zero if all color values are the same
    if color_range == 0:
        normalized_differences = np.zeros_like(color_differences)
    else:
        normalized_differences = color_differences / color_range * 2

    # Step 4: Compute average normalized color difference
    avg_normalized_difference = np.mean(normalized_differences)

    # print("Average Normalized Color Difference (by range):", avg_normalized_difference)
    return avg_normalized_difference


def calc_svd(folder, parameters, use_Vh=None, note="", save_basis=False):

    all_feature, all_feature_name, all_feature_tex, all_Sq, qB = get_all_feature_Sq_data(folder, parameters)

    print("all_feature shape:", np.array(all_feature).shape)
    print("all_Sq shape:", np.array(all_Sq).shape)
    svd = np.linalg.svd(all_Sq)
    print(svd.S)
    print("np.array(svd.U).shape", np.array(svd.U).shape)
    print("np.array(svd.S).shape", np.array(svd.S).shape)
    print("np.array(svd.Vh).shape", np.array(svd.Vh).shape)
    # print(np.linalg.svd(all_Delta_Sq))

    plt.figure(figsize=(6, 6))
    # Subplot for svd.S
    ax00 = plt.subplot(2, 1, 1)
    ax00.plot(range(len(svd.S)), svd.S, "o--", markerfacecolor="none", label="svd.S")
    ax00.set_title("Singular Values (svd.S)")

    # Subplot for svd.U
    ax01 = plt.subplot(2, 1, 2)

    print("np.minimum(svd.Vh[0]), np.maximum(svd.Vh[0])", svd.Vh[0].min(), svd.Vh[0].max())
    ax01.semilogx(qB, svd.Vh[0], "-", markerfacecolor="none", label="svd.Vh[0]")
    ax01.semilogx(qB, svd.Vh[1], "-", markerfacecolor="none", label="svd.Vh[1]")
    ax01.semilogx(qB, svd.Vh[2], "-", markerfacecolor="none", label="svd.Vh[2]")
    ax01.set_title("Left Singular Vectors (svd.Vh[0])")

    plt.tight_layout()
    plt.savefig(f"{folder}/svd{note}.png", dpi=300)
    # plt.show()
    plt.close()

    # save these analyzed data for further easy plotting
    # svd data
    #VhT = np.transpose(svd.Vh)
    if save_basis:
        data = np.column_stack((qB.flatten(), svd.S, svd.Vh[0], svd.Vh[1], svd.Vh[2]))
        column_names = ["qB", "svd.S", "svd.Vh[0]", "svd.Vh[1]", "svd.Vh[2]"]
        np.savetxt(f"{folder}/data_svd{note}.txt", data, delimiter=",", header=",".join(column_names), comments="")

    # Check if svd.Vh is symmetric
    #is_symmetric = np.allclose(svd.Vh, svd.Vh.T, atol=1e-8)
    #print("Is svd.Vh symmetric?", is_symmetric)
    #print("svd.Vh[0],svd.Vh.T[0]",svd.Vh[0],svd.Vh.T[0])
    if use_Vh is not None:
        SqV = np.dot(all_Sq, np.transpose(use_Vh))
        #SqV = np.inner(all_Sq, use_Vh)
    else:
        SqV = np.dot(all_Sq, np.transpose(svd.Vh))
        #SqV = np.inner(all_Sq, svd.Vh)
    print("SqV shape:", SqV.shape)
    print("len(SqV[:,0])", len(SqV[:,0]))
    print("svd.Vh shape", np.array(svd.Vh).shape)
    print("all_Sq shape", np.array(all_Sq).shape)
    print("np.inner(all_Sq[0],svdVh[0])", np.inner(all_Sq[0], svd.Vh[0]))
    print("np.inner(all_Sq[0], np.transpose(svdVh[0]))", np.inner(all_Sq[0], np.transpose(svd.Vh)[0]))
    print("np.dot(all_Sq, np.transpose(svdVh))[0,0]", np.dot(all_Sq, np.transpose(svd.Vh)[0,0]))
    #print("all_Sq[0]", all_Sq[0])
    #print("SqV[0]", SqV[0])
    plt.figure()
    fig = plt.figure(figsize=(2 * len(all_feature_name), 8))
    axs = [fig.add_subplot(2, len(all_feature_name) // 2, i + 1, projection="3d") for i in range(len(all_feature_name))]
    for i in range(len(all_feature_name)):
        #scatter = axs[i].scatter(SqV[:, 0], SqV[:, 1], SqV[:, 2], c=all_feature[:, i], cmap="jet_r", s=1)

        scatter = axs[i].scatter(SqV[:, 0], SqV[:, 1], SqV[:, 2], c=all_feature[:, i], cmap="jet_r", s=1)
        NND = calc_nearest_neighbor_distance(SqV, all_feature[:, i])

        axs[i].set_xlabel("V[0]")
        axs[i].set_ylabel("V[1]")
        axs[i].set_zlabel("V[2]")
        axs[i].set_title(all_feature_tex[i] + f" NND={NND:.2f}")
        axs[i].set_box_aspect([1, 1, 1])  # Set the aspect ratio of the plot
        # Set the same range for each axis
        max_range = np.array([SqV[:, 0].max() - SqV[:, 0].min(), SqV[:, 1].max() - SqV[:, 1].min(), SqV[:, 2].max() - SqV[:, 2].min()]).max() / 2.0
        mid_x = (SqV[:, 0].max() + SqV[:, 0].min()) * 0.5
        mid_y = (SqV[:, 1].max() + SqV[:, 1].min()) * 0.5
        mid_z = (SqV[:, 2].max() + SqV[:, 2].min()) * 0.5
        axs[i].set_xlim(mid_x - max_range, mid_x + max_range)
        axs[i].set_ylim(mid_y - max_range, mid_y + max_range)
        axs[i].set_zlim(mid_z - max_range, mid_z + max_range)
        fig.colorbar(scatter, ax=axs[i], fraction=0.02)
        axs[i].view_init(elev=10.0, azim=-30)

    plt.tight_layout()
    plt.savefig(f"{folder}/svd_projection_scatter_plot{note}.png", dpi=300)
    plt.show()
    plt.close()

    #  svd projection data
    # save svd projection data
    data = np.column_stack((all_feature, SqV[:, 0], SqV[:, 1], SqV[:, 2]))
    column_names = all_feature_name + ["sqv[0]", "sqv[1]", "sqv[2]"]
    np.savetxt(f"{folder}/data_svd_projection{note}.txt", data, delimiter=",", header=",".join(column_names), comments="")

    return svd.Vh


def NND_analysis(folder, invKs, parameter_slices, use_Vh):
    NNDs = []
    for i in range(len(invKs)):
        invK = invKs[i]
        parameters = parameter_slices[i]
        all_feature, all_feature_name, all_feature_tex, all_Sq, qB = get_all_feature_Sq_data(folder, parameters)
        SqV = np.inner(all_Sq, np.transpose(use_Vh))
        NND = []
        for i in range(len(all_feature_name)):
            NND.append(calc_nearest_neighbor_distance(SqV, all_feature[:, i]))
        NNDs.append(NND)

    NNDs = np.array(NNDs)
    plt.figure()
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    for i in range(len(all_feature_name)):
        ax.plot(invKs, NNDs[:, i], "-o", label=all_feature_tex[i])
    ax.set_xlabel(r"$1/K$")
    ax.set_ylabel("NND")
    ax.set_title("Nearest Neighbor Distance")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{folder}/NND_analysis.png", dpi=300)
    plt.show()
    plt.close()
    # save NND data
    data = np.column_stack((invKs, NNDs))
    column_names = ["invK", *all_feature_name]
    np.savetxt(f"{folder}/data_NND_analysis.txt", data, delimiter=",", header=",".join(column_names), comments="")




def calc_Sq_pair_distance_distribution(all_Delta_Sq, max_z, bin_num):
    all_z = np.linspace(0, max_z, bin_num)
    all_Delta_Sq_dis = np.zeros(bin_num)
    all_Delta_Sq_dis[0] = 1 / len(all_Delta_Sq)  # for self distance

    for i in range(len(all_Delta_Sq) - 1):
        for j in range(i + 1, len(all_Delta_Sq)):
            Delta_Sq_dis = np.sqrt(np.sum(np.square(all_Delta_Sq[i] - all_Delta_Sq[j])))
            bin_index = int(Delta_Sq_dis / (max_z / bin_num))
            if bin_index >= bin_num:
                raise ValueError(f"bin_index >= bin_num, z={bin_index/bin_num*max_z}")
            all_Delta_Sq_dis[bin_index] += 2.0 / (len(all_Delta_Sq)) ** 2 / (max_z / bin_num)  # 2.0 for i,j and j,i symmetry, normalize to 1

    return all_Delta_Sq_dis, all_z


def calc_Sq_autocorrelation(mu, all_Delta_Sq, max_z, bin_num):
    # measure the autocorrelation of mu(Delta_Sq)
    all_z = np.linspace(0, max_z, bin_num)
    avg_mu = np.mean(mu)
    avg_mu2 = np.mean(np.square(mu))
    print("np.shape(mu)", np.shape(mu))
    print("np.shape(all_Delta_Sq)", np.shape(all_Delta_Sq))

    print("avg_mu:", avg_mu)
    print("avg_mu2:", avg_mu2)
    print("avg_mu2-avg_mu**2:", avg_mu2 - avg_mu**2)

    avg_mumuz = np.zeros(bin_num)
    avg_mu2z = np.zeros(bin_num)
    avg_muz = np.zeros(bin_num)

    # avg_mumuz[0] = avg_mu2  # for self distance
    # for i in range(len(all_Delta_Sq)-1):
    #    for j in range(i+1, len(all_Delta_Sq)):
    bin_count = np.zeros(bin_num)
    for i in range(len(mu)):
        for j in range(len(mu)):
            Delta_Sq_dis = np.sqrt(np.sum(np.square(all_Delta_Sq[i] - all_Delta_Sq[j])))
            bin_index = int(Delta_Sq_dis / (max_z / bin_num))
            if bin_index >= bin_num:
                raise ValueError(f"bin_index >= bin_num, z={bin_index/bin_num*max_z}")
            avg_muz[bin_index] += mu[i]
            avg_mu2z[bin_index] += mu[i] * mu[i]
            avg_mumuz[bin_index] += mu[i] * mu[j]
            bin_count[bin_index] += 1
    for i in range(bin_num):
        avg_muz[i] /= bin_count[i]
        avg_mu2z[i] /= bin_count[i]
        avg_mumuz[i] /= bin_count[i]

    ac_mu = np.ones(bin_num)
    for i in range(0, bin_num):
        if avg_mu2 - avg_mu**2 == 0:
            ac_mu[i] = 1
        else:
            ac_mu[i] = (avg_mumuz[i] - avg_muz[i] ** 2) / (avg_mu2z[i] - avg_muz[i] ** 2)

    if ac_mu[0] != 1:
        print("ac_mu[0]!=1: ")
        print("avg_mumuz[0]-avg_muz[0]**2,", avg_mumuz[0] - avg_muz[0] ** 2)
        print("(avg_mu2z[0]-avg_muz[0]**2)", (avg_mu2z[0] - avg_muz[0] ** 2))
        print(bin_count)
    ac_mu[0] = 1
    print("ac_mu", ac_mu)
    return ac_mu, all_z


def plot_pddf_acf(folder, parameters, max_z=2, n_bin=100):

    all_feature, all_feature_names, all_Sq2D_flatten, qB = get_all_feature_Sq2D_data(folder, parameters)

    p_z, z = calc_Sq_pair_distance_distribution(all_Sq2D_flatten, max_z, n_bin)

    plt.figure(figsize=(8, 6))
    plt.plot(z, p_z / np.max(p_z), label="p_z/max(p_z)")

    acf_data = []
    for i in range(len(all_feature_names)):
        # pass
        acf_mu, z = calc_Sq_autocorrelation(all_feature[:, i], all_Sq2D_flatten, max_z, n_bin)
        plt.plot(z, acf_mu, label=f"acf_{all_feature_names[i]}")
        acf_data.append(acf_mu)

    plt.xlabel("z")
    plt.ylabel("Value")
    plt.title("Pair Distance Distribution and Autocorrelation")
    plt.legend()
    plt.savefig(f"{folder}/pddf_acf.png", dpi=300)
    plt.close()

    # save these data to file for futher easy plotting

    data = np.column_stack((z, p_z, *acf_data))
    column_names = ["z", "p_z", *["acf_" + feature_name for feature_name in all_feature_names]]
    np.savetxt(f"{folder}/data_pddf_acf.txt", data, delimiter=",", header=",".join(column_names), comments="")


def GaussianProcess_optimization(folder, parameters_train, feature_to_train, note=""):
    all_feature, all_feature_name, all_feature_tex, all_Sq, qB = get_all_feature_Sq_data(folder, parameters_train)
    grid_size = 30

    theta_per_feature = {
        # "L": (np.logspace(-2, 1, grid_size), np.logspace(-5, -2, grid_size)),
        "kappa": (np.logspace(-2, 0, grid_size), np.logspace(-4, -2, grid_size)),
        "A": (np.logspace(-2, 0, grid_size), np.logspace(-2, 0 , grid_size)),
        #"invK": (np.logspace(-2, 0, grid_size), np.logspace(-3, -1, grid_size)),
        "R2": (np.logspace(-1, 1, grid_size), np.logspace(-5, -3, grid_size)),
        "Rg2": (np.logspace(0, 1, grid_size), np.logspace(-11, -10, grid_size)),
    }

    # feature normalization
    all_feature_mean = np.mean(all_feature, axis=0)
    all_feature_std = np.std(all_feature, axis=0)
    all_feature = (all_feature - all_feature_mean) / all_feature_std
    all_gp_per_feature = {}
    plt.figure()
    fig, axs = plt.subplots(1, len(all_feature_name), figsize=(6 * len(all_feature_name), 6))
    for feature_name, (theta0, theta1) in theta_per_feature.items():
        if feature_name not in feature_to_train:
            continue
        print("training: ", feature_name)
        feature_index = all_feature_name.index(feature_name)

        F_learn = all_Sq

        # witout theta optimization
        kernel = RBF(1) + WhiteKernel(1)
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, optimizer=None).fit(F_learn, all_feature[:, feature_index])
        # print(" all_feature[:, feature_index]", all_feature[:, feature_index])

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        # calc Log likelihood
        ax = axs[all_feature_name.index(feature_name)]
        Theta0, Theta1 = np.meshgrid(theta0, theta1)
        LML = [[0 for j in range(Theta0.shape[1])] for i in range(Theta0.shape[0])]
        for i in range(Theta0.shape[0]):
            for j in range(Theta0.shape[1]):
                LML[i][j] = gp.log_marginal_likelihood(np.log([Theta0[i, j], Theta1[i, j]]))
                print(f"Calculating LML: i={i}/{Theta0.shape[0]}, j={j}/{Theta0.shape[1]}, LML={LML[i][j]}", end="\r")
        # reason for np.log here is the theta is log-transformed hyperparameters (https://github.com/scikit-learn/scikit-learn/blob/5491dc695/sklearn/gaussian_process/kernels.py#L1531) line (289)

        ax.contour(Theta0, Theta1, LML, levels=200)
        # find optimized theta0, theta1, using the above contour as guidanve
        kernel = RBF(theta0[grid_size // 2], (theta0[0], theta0[-1])) + WhiteKernel(theta1[grid_size // 2], (theta1[0], theta1[-1]))
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer=10).fit(F_learn, all_feature[:, feature_index])
        all_gp_per_feature[feature_name] = gp

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        ax.plot([gp_theta[0]], [gp_theta[1]], "x", color="red", markersize=10, markeredgewidth=2, label=r"l=%.2e, $\sigma$=%.2e" % (gp_theta[0], gp_theta[1]))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"theta0: l")
        ax.set_ylabel(r"theta1: $\sigma$")
        feature_name_legend = feature_name
        ax.set_title(f"Log Marginal Likelihood for {feature_name_legend}")
        ax.legend()

        data = np.column_stack(([gp_theta[0]] * len(theta0), [gp_theta[1]] * len(theta1), theta0, theta1, np.array(LML).T))
        column_names = ["gp_theta0", "gp_theta1", "theta0", "theta1", "LML"]
        np.savetxt(f"{folder}/data_{feature_name}_LML{note}.txt", data, delimiter=",", header=",".join(column_names), comments="")
        with open(f"{folder}/gp_{feature_name}.pkl", "wb") as f:
            pickle.dump(gp, f)

    # Save average and standard deviation per feature
    avg_std_data = np.column_stack((all_feature_name, all_feature_mean, all_feature_std))
    column_names = ["Feature", "Mean", "Std"]
    np.savetxt(f"{folder}/data_feature_avg_std{note}.txt", avg_std_data, delimiter=",", header=",".join(column_names), comments="", fmt="%s")

    plt.tight_layout()
    plt.savefig(f"{folder}/LML_subplots{note}.png", dpi=300)
    # plt.show()
    plt.close()

    # return trained GPR
    return all_feature_mean, all_feature_std, all_gp_per_feature


def read_gp_and_feature_stats(folder, note=""):
    all_feature_names = ["L", "kappa", "A", "invK", "R2", "Rg2"]
    all_feature_mean = np.genfromtxt(f"{folder}/data_feature_avg_std{note}.txt", delimiter=",", skip_header=1, usecols=1)
    all_feature_std = np.genfromtxt(f"{folder}/data_feature_avg_std{note}.txt", delimiter=",", skip_header=1, usecols=2)
    all_gp_per_feature = {}
    for feature_name in all_feature_names:
        if os.path.exists(f"{folder}/gp_{feature_name}.pkl"):
            with open(f"{folder}/gp_{feature_name}.pkl", "rb") as f:
                all_gp_per_feature[feature_name] = pickle.load(f)
    return all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature


def GaussianProcess_prediction(folder, parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature, feature_to_test, note=""):
    all_feature, all_feature_name, all_feature_tex, all_Sq, qB = get_all_feature_Sq_data(folder, parameters_test)

    plt.figure()

    fig, axs = plt.subplots(1, len(all_feature_name), figsize=(6 * len(all_feature_name), 6))
    for feature_name, gp in all_gp_per_feature.items():
        if feature_name not in feature_to_test:
            continue
        feature_index = all_feature_name.index(feature_name)
        Y = all_feature[:, feature_index]

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)  # gp.kernel_.theta return log transformed theta
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        Y_predict, Y_predict_err = gp.predict(all_Sq, return_std=True)
        # print("np.shape(test_data[:, 0])", np.shape(test_data[:, 0]))
        print("np.shape(all_Sq2D_flatten)", np.shape(all_Sq))
        print("np.shape(Y_predict)", np.shape(Y_predict))

        Y_predict = Y_predict * all_feature_std[feature_index] + all_feature_mean[feature_index]
        Y_predict_err = Y_predict_err * all_feature_std[feature_index]

        axs[feature_index].errorbar(Y, Y_predict, yerr=Y_predict_err, marker="o", markerfacecolor="none", markersize=3, linestyle="none")
        axs[feature_index].plot(Y, Y, "--")
        min_val = min(np.min(Y), np.min(Y_predict - Y_predict_err))
        max_val = max(np.max(Y), np.max(Y_predict + Y_predict_err))
        axs[feature_index].set_xlim(min_val, max_val)
        axs[feature_index].set_ylim(min_val, max_val)
        axs[feature_index].set_xlabel(f"{all_feature_tex[feature_index]}")
        axs[feature_index].set_ylabel(f"{all_feature_tex[feature_index]} Prediction")

        # save data to file
        data = np.column_stack((Y, Y_predict, Y_predict_err))
        column_names = [feature_name, "ML predicted", "ML predicted uncertainty"]
        np.savetxt(f"{folder}/data_{feature_name}_prediction{note}.txt", data, delimiter=",", header=",".join(column_names), comments="")

    plt.savefig(f"{folder}/prediction{note}.png", dpi=300)
    plt.close()


def ax_fit(x, a):
    return a * x


def fit_Rg2(q, Sq):
    popt, pcov = curve_fit(ax_fit, q**2 / 3, (1 - Sq))
    perr = np.sqrt(np.diag(pcov))
    return popt[0], perr[0]


def calc_Sq_fitted_Rg2(folder, parameters_test, all_feature_names):
    segment_type, all_feature, all_feature_names, all_Sq, all_Sq_err, q = read_Sq_data(folder, parameters_test)

    MC_Rg2 = all_feature[:, all_feature_names.index("Rg2")]
    # qfns = [10,20,30,40]
    qfns = [50, 55, 60, 65, 70]
    Rg2s = []
    Rg2_errs = []
    plt.figure()
    for qfn in qfns:
        Rg2s.append([])
        Rg2_errs.append([])
        for i in range(len(all_Sq)):
            Rg2, Rg2_err = fit_Rg2(q[:qfn], all_Sq[i][:qfn])
            Rg2s[-1].append(Rg2)
            Rg2_errs[-1].append(Rg2_err)

        plt.scatter(MC_Rg2, Rg2s[-1], alpha=0.5, label=f"qf={q[qfn-1]}")
    plt.plot(MC_Rg2, MC_Rg2, "k--")
    plt.xlabel("MC Rg2")
    plt.ylabel("Fitted Rg2")
    plt.legend()
    plt.savefig(f"{folder}/{segment_type}_Rg2_fit.png", dpi=300)
    plt.close()

    data = np.column_stack(([MC_Rg2] + Rg2s))
    column_names = ["MC Rg2", "fitted Rg2"]
    np.savetxt(f"{folder}/data_{segment_type}_fitted_Rg2.txt", data, delimiter=",", header=",".join(column_names), comments="")
