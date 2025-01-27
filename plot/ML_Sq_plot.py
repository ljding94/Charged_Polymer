import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib import rc
import os
from conformation_plot import get_obs_data
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def calc_Sq_discrete_infinite_thin_rod(q, L):
    # numereical calculation
    Sq = [1.0 / L for i in range(len(q))]
    for k in range(len(q)):
        Sqk = 0
        qk = q[k]
        for i in range(L - 1):
            for j in range(i + 1, L):
                Sqk += 2.0 * np.sin(qk * (i - j)) / (qk * (i - j)) / (L * L)
        Sq[k] += Sqk
    return np.array(Sq)


def plot_Sq_example(tex_lw=240.71031, ppi=72):
    print("\nplotting Sq figure\n")

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222, sharex=ax0)
    ax2 = fig.add_subplot(223, sharex=ax0)
    ax3 = fig.add_subplot(224, sharex=ax0)

    folder = "../data/20250108"

    # plot various kappa
    L = 500
    kappa = 30.0
    A = 5.0
    invK = 3.0
    Sq_rod = None
    for kappa in [10, 30, 50]:
        data = get_obs_data(folder, L, kappa, A, invK)
        Sq = data[0, 10:]
        q = data[2, 10:]
        if Sq_rod is None:
            Sq_rod = calc_Sq_discrete_infinite_thin_rod(q, L)
            ax0.loglog(q, Sq_rod, "k--", lw=1, label=f"rod")
        ax0.loglog(q, Sq, "-", lw=1, label=f"{kappa}")
        ax1.semilogx(q, Sq / Sq_rod, "-", lw=1, label=f"{kappa}")

    kappa = 30.0

    # ax0.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax0.set_ylabel(r"$S(q)$", fontsize=9, labelpad=0)
    ax0.tick_params(which="both", direction="in", top="on", right="on", labelbottom=None, labelsize=7)
    ax0.legend(title=r"$\kappa$", loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # ax1.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax1.set_ylabel(r"$S(q)/S_{rod}(q)$", fontsize=9, labelpad=0)
    ax1.tick_params(which="both", direction="in", top="on", right="on", labelbottom=None, labelsize=7)
    ax1.legend(title=r"$\kappa$", loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax1.yaxis.set_major_locator(MultipleLocator(0.4))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.2))

    # plot various A
    for A in [1, 3, 5, 7]:
        data = get_obs_data(folder, L, kappa, A, invK)
        Sq = data[0, 10:]
        q = data[2, 10:]
        ax2.semilogx(q, Sq / Sq_rod, "-", lw=1, label=f"{A}")
    A = 5.0

    ax2.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax2.set_ylabel(r"$S(q)/S_{rod}(q)$", fontsize=9, labelpad=0)
    ax2.tick_params(which="both", direction="in", top="on", right="on", labelsize=7)
    ax2.legend(title=r"$A$", loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax2.yaxis.set_major_locator(MultipleLocator(0.4))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.2))

    # plot various invK
    for invK in [1, 3, 5, 7]:
        data = get_obs_data(folder, L, kappa, A, invK)
        Sq = data[0, 10:]
        q = data[2, 10:]
        ax3.semilogx(q, Sq / Sq_rod, "-", lw=1, label=f"{invK}")
    invK = 3.0
    ax3.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax3.set_ylabel(r"$S(q)/S_{rod}(q)$", fontsize=9, labelpad=0)
    ax3.tick_params(which="both", direction="in", top="on", right="on", labelsize=7)
    ax3.legend(title=r"$\lambda_D$", loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax3.yaxis.set_major_locator(MultipleLocator(0.4))
    ax3.yaxis.set_minor_locator(MultipleLocator(0.2))

    ax0.text(0.7, 0.15, r"$(a)$", fontsize=9, transform=ax0.transAxes, color="black")
    ax1.text(0.7, 0.15, r"$(b)$", fontsize=9, transform=ax1.transAxes, color="black")
    ax2.text(0.7, 0.15, r"$(c)$", fontsize=9, transform=ax2.transAxes, color="black")
    ax3.text(0.7, 0.15, r"$(d)$", fontsize=9, transform=ax3.transAxes, color="black")

    plt.tight_layout(pad=0.2)
    plt.savefig("./figures/Sq_example.pdf", format="pdf")
    plt.savefig("./figures/Sq_example.png", dpi=300)
    plt.show()
    plt.close()


def get_all_feature_Sq_data(folder, nrun=15000, filter_invK=None):

    all_L, all_kappa, all_A, all_invK = [], [], [], []
    all_R2, all_Rg2 = [], []  # R2, Rg2 related
    all_Sq = []
    qB = []
    for i in range(nrun):
        filename = f"{folder}/obs_random_run{i}.csv"
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
        all_R2.append(R2)
        all_Rg2.append(Rg2 / L)

        Sq = data[0, 10:]
        all_Sq.append(Sq)
        qB = data[2, 10:]

    all_feature = np.array([all_L, all_kappa, all_A, all_invK, all_R2, all_Rg2])
    all_feature_name = ["L", "kappa", "A", "invK", "R2", "Rg2"]
    all_feature_tex = [r"L", r"$\kappa$", "A", r"$1/K$", r"$R^2$", r"$R_g^2/L$"]
    all_Sq = np.array(all_Sq)
    all_Sq = np.log(all_Sq)
    qB = np.array(qB)
    return all_feature.T, all_feature_name, all_feature_tex, all_Sq, qB


def plot_SVD_basis(tex_lw=240.71031, ppi=72):
    folder = "../data/20250114_rand"

    data = np.loadtxt(f"{folder}/data_svd_all.txt", skiprows=1, delimiter=",", unpack=True)
    q, S, V0, V1, V2 = data
    V0, V1, V2 = -V0, -V1, -V2

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    ax1.semilogx(range(1, len(S) + 1), S, "x", color="royalblue", mfc="None")
    ax1.semilogx(range(1, 4), S[:3], "o", color="red", mfc="None")
    ax1.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax1.set_xlabel(r"SVR", fontsize=9, labelpad=0)
    ax1.set_ylabel(r"$\Sigma$", fontsize=9, labelpad=0)
    ax1.yaxis.set_major_locator(plt.MultipleLocator(500))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(250))

    # ax2 = ax1.inset_axes([0.2, 0.2, 0.9, 0.9])

    ax2.semilogx(q, V0, lw=1, label=r"$V_0$")
    ax2.semilogx(q, V1, lw=1, label=r"$V_1$")
    ax2.semilogx(q, V2, lw=1, label=r"$V_2$")
    ax2.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax2.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax2.set_ylabel(r"$V$", fontsize=9, labelpad=-5)
    # ax2.legend(loc="center right", bbox_to_anchor=(1.0,0.65),ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)
    ax2.legend(loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)
    ax2.set_xlim(4e-4, 1.1)

    data = get_obs_data("../data/20250108", 500, 10, 5, 3)
    Sq = data[0, 10:]
    q = data[2, 10:]

    print("lnSq0/V0", np.inner(np.log(Sq), V0))
    print("lnSq1/V1", np.inner(np.log(Sq), V1))
    print("lnSq2/V2", np.inner(np.log(Sq), V2))

    lnSq0 = np.inner(np.log(Sq), V0) * V0
    lnSq1 = np.inner(np.log(Sq), V1) * V1
    lnSq2 = np.inner(np.log(Sq), V2) * V2
    # ax3.loglog(q, Sq, linestyle=(0, (2, 2)), color="red", lw=1.5, label=r"$S(QB)$")
    ax3.loglog(q, Sq, ls="-", color="red", lw=1.5, label=r"$S(q)$")
    ax3.loglog(q, np.exp(lnSq0 + lnSq1 + lnSq2), linestyle=(0, (2, 2)), lw=1.5, color="blue", label=r"$S_0\odot S_1\odot S_2$")
    ax3.loglog(q, np.exp(lnSq0), "-", lw=0.75, label=r"$S_0$")
    ax3.loglog(q, np.exp(lnSq1), "-", lw=0.75, label=r"$S_1$")
    ax3.loglog(q, np.exp(lnSq2), "-", lw=0.75, label=r"$S_2$")

    ax3.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax3.set_xlabel(r"$q$", fontsize=9, labelpad=0)
    ax3.set_ylabel(r"$S(q)$", fontsize=9, labelpad=0)
    ax3.legend(loc="lower left", ncol=2, columnspacing=0.5, handlelength=1.5, handletextpad=0.5, frameon=False, fontsize=9)

    ax1.text(0.8, 0.2, r"$(a)$", fontsize=9, transform=ax1.transAxes, color="black")
    ax2.text(0.8, 0.2, r"$(b)$", fontsize=9, transform=ax2.transAxes, color="black")
    ax3.text(0.9, 0.2, r"$(c)$", fontsize=9, transform=ax3.transAxes, color="black")

    plt.tight_layout(pad=0.5)
    plt.savefig("figures/SVD_basis_Sq.pdf", format="pdf")
    plt.savefig("figures/SVD_basis_Sq.png", dpi=300)
    plt.show()
    plt.close()


def plot_SVD_feature(tex_lw=240.71031, ppi=72):
    print("\nplotting SVD feature figure\n")
    folder = "../data/20250114_rand"
    data = np.loadtxt(
        f"{folder}/data_svd_projection_all.txt",
        skiprows=1,
        delimiter=",",
        unpack=True,
    )
    L, kappa, A, invK, R2, Rg2, sqv0, sqv1, sqv2 = data
    sqv0, sqv1, sqv2 = -sqv0, -sqv1, -sqv2

    feature_tex_to_data = {r"$R^2/L^2$": R2, r"$R_g^2/L^2$": Rg2, r"$\kappa$": kappa, r"$A$": A}

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax0 = fig.add_subplot(221, projection="3d")
    ax1 = fig.add_subplot(222, projection="3d")
    ax2 = fig.add_subplot(223, projection="3d")
    ax3 = fig.add_subplot(224, projection="3d")

    cbar_major_locator = [0.4, 0.02, 20, 3]
    axis_major_locator = [1, 2, 1]
    axis_lim = [(21.5, 25), (-1, 4), (-1, 2.5)]
    axs = [ax0, ax1, ax2, ax3]
    for i, (key, data) in enumerate(feature_tex_to_data.items()):
        scatter = axs[i].scatter(sqv0, sqv1, sqv2, s=0.5, c=data, cmap="rainbow", rasterized=True)
        axs[i].set_xlabel(r"$FV_0$", fontsize=9, labelpad=-10, rotation=0)
        axs[i].set_ylabel(r"$FV_1$", fontsize=9, labelpad=-10, rotation=0)
        axs[i].set_zlabel(r"$FV_2$", fontsize=9, labelpad=-10, rotation=0)
        axs[i].tick_params("x", direction="in", which="both", labelsize=7, pad=-5, rotation=0)
        axs[i].tick_params("y", direction="in", which="both", labelsize=7, pad=-5)
        axs[i].tick_params("z", direction="in", which="both", labelsize=7, pad=-5)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[0]))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[1]))
        axs[i].zaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[2]))
        axs[i].set_xlim(axis_lim[0])
        axs[i].set_ylim(axis_lim[1])
        axs[i].set_zlim(axis_lim[2])

        cbar = fig.colorbar(scatter, ax=axs[i], fraction=0.02, pad=-0.04)  # , location="left")
        cbar.ax.tick_params(direction="in", which="both", labelsize=7, pad=1.0)
        cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(cbar_major_locator[i]))
        cbar.ax.yaxis.set_minor_locator(plt.MultipleLocator(cbar_major_locator[i] * 0.5))
        cbar.ax.set_title(key, fontsize=9)
        axs[i].grid(True, which="minor")
        axs[i].view_init(elev=-120, azim=120)

    textlabel = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]
    for i in range(len(axs)):
        axs[i].text2D(0.2, 0.85, textlabel[i], fontsize=9, transform=axs[i].transAxes, color="black")

    plt.tight_layout(pad=0.3, w_pad=0.5, h_pad=1)
    plt.subplots_adjust(left=-0.07, bottom=0.07)
    plt.savefig("figures/SVD_feature.pdf", format="pdf", dpi=300)
    plt.savefig("figures/SVD_feature.png", dpi=300)
    plt.show()
    plt.close()


def plot_SVD_NND(tex_lw=240.71031, ppi=72):
    print("\nplotting SVD NND figure\n")
    folder = "../data/20250114_rand"
    data = np.loadtxt(
        f"{folder}/data_NND_analysis.txt",
        skiprows=1,
        delimiter=",",
        unpack=True,
    )
    invK, NND_L, NND_kappa, NND_A, NND_invK, NND_R2, NND_Rg2 = data

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax00 = plt.subplot2grid((5, 3), (0, 0), projection="3d", rowspan=2)
    ax01 = plt.subplot2grid((5, 3), (0, 1), projection="3d", rowspan=2)
    ax02 = plt.subplot2grid((5, 3), (0, 2), projection="3d", rowspan=2)
    # ax03 = plt.subplot2grid((3, 4), (0, 3), projection="3d")
    ax10 = plt.subplot2grid((5, 3), (2, 0), colspan=3, rowspan=3)

    ax10.plot(invK, NND_R2, "-^", lw=1, mfc="none", ms=4, label=r"$R^2/L^2$")
    ax10.plot(invK, NND_Rg2, "-v", lw=1, mfc="none", ms=4, label=r"$R_g^2/L^2$")
    ax10.plot(invK, NND_A, "-o", lw=1, mfc="none", ms=4, label=r"$A$")
    ax10.plot(invK, NND_kappa, "-s", lw=1, mfc="none", ms=4, label=r"$\kappa$")

    ax10.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax10.set_ylabel(r"$D_{NN}$", fontsize=9, labelpad=0)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelsize=7)
    ax10.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # get svd per invK
    PinvKs = [1, 3, 5]
    axs = [ax00, ax01, ax02]
    axis_major_locator = [1, 2, 1]
    axis_lim = [(21.5, 25), (-1, 4), (-1, 2.5)]
    for i in range(len(PinvKs)):
        invK = PinvKs[i]
        print(f"invK: {invK}, invK-0.5:{invK-0.5:.1f}, invK+0.5:{invK+0.5:.1f}")
        data = np.loadtxt(f"{folder}/data_svd_projection_invK_{invK-0.5:.1f}_{invK+0.5:.1f}.txt", skiprows=1, delimiter=",", unpack=True)
        L, kappa, A, vinvK, R2, Rg2, sqv0, sqv1, sqv2 = data
        sqv0, sqv1, sqv2 = -sqv0, -sqv1, -sqv2

        scatter = axs[i].scatter(sqv0, sqv1, sqv2, s=0.5, c=A, cmap="rainbow", rasterized=True)
        # axs[i].set_xlabel(r"$FV_0$", fontsize=9, labelpad=-6, rotation=0)
        # axs[i].set_ylabel(r"$FV_1$", fontsize=9, labelpad=-6, rotation=0)
        # axs[i].set_zlabel(r"$FV_2$", fontsize=9, labelpad=-6, rotation=0)
        axs[i].tick_params("x", direction="in", which="both", labelsize=7, pad=-5, rotation=-15)
        axs[i].tick_params("y", direction="in", which="both", labelsize=7, pad=-2)
        axs[i].tick_params("z", direction="in", which="both", labelsize=7, pad=-2)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[0]))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[1]))
        axs[i].zaxis.set_major_locator(plt.MultipleLocator(axis_major_locator[2]))
        axs[i].set_xlim(axis_lim[0])
        axs[i].set_ylim(axis_lim[1])
        axs[i].set_zlim(axis_lim[2])

        # axs[i].set_box_aspect([1, 1, 1.5])
        axs[i].tick_params(axis="both", which="both", length=0, labelbottom=False, labelleft=False)
        axs[i].grid(True, which="minor")
        axs[i].view_init(elev=-120, azim=120)
        axs[i].text2D(0.6, -0.01, rf"$\lambda_D$={invK}", fontsize=9, transform=axs[i].transAxes, color="black")

    cbar = fig.colorbar(scatter, ax=ax02, fraction=0.03, pad=0)
    cbar.ax.tick_params(direction="in", which="both", labelsize=7, pad=0.8)
    cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(3))
    cbar.ax.yaxis.set_minor_locator(plt.MultipleLocator(1.5))
    cbar.ax.set_title(r"$A$", fontsize=9)

    ax00.text2D(0.05, 0.1, r"$(a)$", fontsize=9, transform=ax00.transAxes, color="black")
    ax10.text(0.02, 0.4, r"$(b)$", fontsize=9, transform=ax10.transAxes, color="black")

    plt.tight_layout(pad=0.5)
    plt.savefig("figures/SVD_NND.pdf", format="pdf", dpi=300)
    plt.savefig("figures/SVD_NND.png", dpi=300)
    plt.show()
    plt.close()


def plot_LML(tex_lw=240.71031, ppi=72):
    print("\nplotting LML figure\n")
    folder = "../data/20250114_rand"

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((2, 2), (0, 0))
    ax01 = plt.subplot2grid((2, 2), (0, 1))
    ax10 = plt.subplot2grid((2, 2), (1, 0))
    ax11 = plt.subplot2grid((2, 2), (1, 1))

    features = ["R2", "Rg2", "kappa", "A"]
    feature_texs = [r"$R^2/L^2$", r"$R_g^2/L^2$", r"$\kappa$", r"$A$"]
    axs = [ax00, ax01, ax10, ax11]
    for i in range(len(features)):
        feature = features[i]
        ax = axs[i]
        note = "_A" if feature == "A" else "_all"
        data = np.loadtxt(f"{folder}/data_{feature}_LML{note}.txt", skiprows=1, delimiter=",", unpack=True)
        gp_theta0, gp_theta1, theta0, theta1, LML = data[0][0], data[1][0], data[2], data[3], data[4:]
        print(f"gp_theta0={gp_theta0}, gp_theta1={gp_theta1}")
        Theta0, Theta1 = np.meshgrid(theta0, theta1)

        ax.contour(Theta0, Theta1, LML, linewidths=1, levels=200)
        ax.plot([gp_theta0], [gp_theta1], "x", color="k", markersize=5, markeredgewidth=1)  # , label=r"l=%.2e, $\sigma$=%.2e" % (gp_theta0, gp_theta1))
        ax.set_xlabel(r"$l$", fontsize=9, labelpad=0)
        ax.set_ylabel(r"$\sigma$", fontsize=9, labelpad=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.text(0.2, 0.8, feature_texs[i], fontsize=9, transform=ax.transAxes, color="black")
        ax.tick_params(labelsize=7, which="both", direction="in", top="on", right="on")  # , labelleft=False, labelbottom=False)
        axs[i].minorticks_off()

    ax00.text(0.8, 0.1, r"$(a)$", fontsize=9, transform=ax00.transAxes)
    ax01.text(0.8, 0.1, r"$(b)$", fontsize=9, transform=ax01.transAxes)
    ax10.text(0.8, 0.1, r"$(c)$", fontsize=9, transform=ax10.transAxes)
    ax11.text(0.8, 0.1, r"$(d)$", fontsize=9, transform=ax11.transAxes)


    plt.tight_layout(pad=0.5)
    # plt.subplots_adjust(left=0.2, bottom=0.1)
    plt.savefig("figures/LML.pdf", format="pdf", dpi=300)
    plt.savefig("figures/LML.png", dpi=300)
    plt.show()
    plt.close()


def plot_GPR(tex_lw=240.71031, ppi=72):
    pass
    print("\nplotting GPR  figure\n")
    folder = "../data/20250114_rand"

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)

    features = ["R2", "Rg2", "kappa", "A"]
    feature_texs = [r"$R^2/L^2$", r"$R_g^2/L^2$", r"$\kappa$", r"$A$"]
    major_locator = [0.3, 0.02, 15, 3]
    axs = [ax0, ax1, ax2, ax3]
    for i in range(len(features)):
        feature = features[i]
        ax = axs[i]
        note = "_A" if feature == "A" else "_all"
        data = np.loadtxt(f"{folder}/data_{feature}_prediction{note}.txt", skiprows=1, delimiter=",", unpack=True)
        mu, mu_pred, mu_pred_err = data[0], data[1], data[2]

        axs[i].scatter(mu, mu_pred, s=0.5, marker=".", c="royalblue")
        err = np.average(np.abs(mu - mu_pred) / np.maximum(mu, mu_pred))
        axs[i].text(0.5, 0.3, f"Err={err:.2f}", transform=axs[i].transAxes, fontsize=9)
        axs[i].plot(mu, mu, color="gray", linestyle="--", lw=0.25, alpha=0.5)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        #axs[i].xaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        #axs[i].yaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].grid(True, which="major", linestyle="--", linewidth=0.5)
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        # axs[i].legend(title = feature_tex, fontsize=9, loc="upper left")
        axs[i].text(0.2, 0.6, feature_texs[i], transform=axs[i].transAxes, fontsize=9)

    axall = fig.add_subplot(111, frameon=False)
    axall.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    axall.set_xlabel("MC References", fontsize=9, labelpad=-3)
    axall.set_ylabel("ML Inversion", fontsize=9, labelpad=-3)
    # ax13.set_ylabel("ML Inversion", fontsize=9, labelpad=0)

    ax0.text(0.8, 0.1, r"$(a)$", fontsize=9, transform=ax0.transAxes)
    ax1.text(0.8, 0.1, r"$(b)$", fontsize=9, transform=ax1.transAxes)
    ax2.text(0.8, 0.1, r"$(c)$", fontsize=9, transform=ax2.transAxes)
    ax3.text(0.8, 0.1, r"$(d)$", fontsize=9, transform=ax3.transAxes)


    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(left=0.125, bottom=0.125)
    plt.savefig("figures/GPR_prediction.pdf", format="pdf", dpi=300)
    plt.savefig("figures/GPR_prediction.png", dpi=300)
    plt.show()
    plt.close()