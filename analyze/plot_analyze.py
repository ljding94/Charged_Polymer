import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib import colormaps as cm
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
import os


def plot_polymer_config(filename, finfo, show=False):
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    # ax.plot(x, y, z, "-o", color="gray", mfc="black", mec="black")
    ax.plot(x, y, z, "-", color="gray")
    ax.plot([x[0]], [y[0]], [z[0]], "o", color="red", label="start")
    ax.plot([x[-1]], [y[-1]], [z[-1]], "o", color="blue", label="end")

    # ymax = max(np.max(y), -np.min(y))
    # zmax = max(np.max(z), -np.min(z))
    ax.set_xlim(np.min(x) - 1, np.max(x) + 1)
    ax.set_ylim(np.min(y) - 1, np.max(y) + 1)
    ax.set_zlim(np.min(z) - 1, np.max(z) + 1)

    # ax.set_ylim(-ymax-1, ymax+1)
    # ax.set_zlim(-zmax-1, zmax+1)

    ax.legend()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(finfo)
    """
    cmap = cm.get_cmap("jet_r")
    norm = Normalize(vmin=0, vmax=0.5*np.pi)
    deg = np.arccos(np.abs(n[:, 2]))
    for i in range(len(r)):
        # ax.plot(np.append(beadsrx[1:], beadsrx[1]), np.append(beadsry[1:], beadsry[1]), np.append(beadsrz[1:], beadsrz[1]), "-", color="gray")
        pass
    cbar = plt.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")), ax=ax, ticks=[0, np.pi/6, np.pi/3, np.pi/2])
    cbar.ax.set_yticklabels([r"$0$", r"$\pi/6$", r"$\pi/3$", r"$\pi/2$"])
    cbar.ax.tick_params(direction="in")
    cbar.ax.set_title(r"$\arccos{|\hat{\mathrm{n}}\cdot\hat{\mathrm{z}}|}$")

    cbar.set_label('Z')
    """
    """
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.zaxis.set_major_locator(plt.MultipleLocator(1))
    ax.zaxis.set_minor_locator(plt.MultipleLocator(0.5))
    """
    ax.set_aspect("equal")
    plt.savefig(filename.replace(".csv", ".png"))
    if show:
        plt.show()
    plt.close()


def plot_multi_config(folder, finfo):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for i in range(30):
        filename = f"{folder}/{finfo}/config_{i}.csv"
        if not os.path.exists(filename):
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        # ax.plot(data[:, 0], data[:, 1], data[:, 2], "o-")
        ax.plot(data[:, 0], data[:, 1], data[:, 2], "-")
    # ax.set_xlim(-15, 15)
    # ax.set_ylim(-15, 15)
    # ax.set_zlim(-15, 15)
    ax.set_aspect("equal")
    plt.savefig(f"{folder}/config_{finfo}.png")
    plt.close()


def plot_MC_step(filename, finfo, show=False):
    data = np.genfromtxt(filename, delimiter=",", skip_header=4)
    E, Tb, X, Y, Z, R = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]

    fig, axs = plt.subplots(6, 1, figsize=(5, 15))

    axs[0].plot(range(len(Tb)), Tb, "*")
    axs[0].set_xlabel("MC sweep (1000L step per sweep)")
    axs[0].set_ylabel("Tb")

    axs[1].plot(range(len(X)), X, "*")
    axs[1].set_xlabel("MC sweep (1000L step per sweep)")
    axs[1].set_ylabel("X")

    axs[2].plot(range(len(Y)), Y, "*")
    axs[2].set_xlabel("MC sweep (1000L step per sweep)")
    axs[2].set_ylabel("Y")

    axs[3].plot(range(len(Z)), Z, "*")
    axs[3].set_xlabel("MC sweep (1000L step per sweep)")
    axs[3].set_ylabel("Z")

    axs[4].plot(range(len(R)), R, "*")
    axs[4].set_xlabel("MC sweep (1000L step per sweep)")
    axs[4].set_ylabel("R")

    axs[5].plot(range(len(E)), E, "*")
    axs[5].set_xlabel("MC sweep (1000L step per sweep)")
    axs[5].set_ylabel("Energy")
    plt.tight_layout()
    plt.savefig(filename.replace(".csv", ".png"))
    if show:
        plt.show()
    plt.close()


def plot_obs(folder, finfos, parameters, xparam, tag=""):
    fig, axs = plt.subplots(1, 5, figsize=(15, 3))
    xpar = []
    all_R2, all_R2_err = [], []
    all_Rg2, all_Rg2_err = [], []
    all_Sq = []
    all_qB = []
    all_tts = []
    all_spB = []

    for i in range(len(finfos)):
        finfo = finfos[i]
        L, kappa, A, invK = parameters[i]
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        E, E_B, E_HSY, R2, Rg2 = data[0, 5], data[0, 6], data[0, 7], data[0, 8], data[0, 9]
        E_err, E_B_err, E_HSY_err, R2_err, Rg2_err = data[1, 5], data[1, 6], data[1, 7], data[1, 8], data[1, 9]

        Sq = data[0, 10:]
        qB = data[2, 10:]
        tts = data[3, 10:]
        spB = data[5, 10:]

        if xparam == "L":
            xpar.append(L)
        elif xparam == "kappa":
            xpar.append(kappa)
        elif xparam == "A":
            xpar.append(A)
        elif xparam == "invK":
            xpar.append(invK)
        else:
            raise ValueError(f"Unknown xparam={xparam}")

        all_R2.append(R2/L)
        all_R2_err.append(R2_err/L)
        all_Rg2.append(Rg2/L)
        all_Rg2_err.append(Rg2_err/L)

        all_Sq.append(Sq)
        all_qB.append(qB)
        all_tts.append(tts)
        all_spB.append(spB)

    axs[0].errorbar(xpar, all_R2, yerr=all_R2_err, marker="+", label=r"$R^2/L$")
    axs[1].errorbar(xpar, all_Rg2, yerr=all_Rg2_err, marker="x", label=r"$R_g^2/L$")

    axs[0].set_xlabel(xparam)
    axs[0].set_ylabel("")
    axs[0].legend()

    axs[1].set_xlabel(xparam)
    axs[1].set_ylabel("")
    axs[1].legend()

    Sq_rod = calc_Sq_discrete_infinite_thin_rod(all_qB[0], parameters[0][0])
    for i in range(len(all_Sq)):
        axs[2].loglog(all_qB[i], all_Sq[i], "-", label=f"{xpar[i]}")
    axs[2].loglog(all_qB[0], Sq_rod, "k--", label="rod")
    axs[2].set_xlabel("qB")
    axs[2].set_ylabel("S(qB)")
    axs[2].legend(title=xparam, ncol=2)

    all_alpha, all_alpha_err = [], []
    all_lam1, all_lam1_err = [], []
    all_lam2, all_lam2_err = [], []
    all_lam0 = []
    for i in range(len(all_tts)):
        #if np.isnan(all_tts[i]).any():
        #    print(f"NaN values found in all_tts[{i}] with parameters: {parameters[i]}")
        L, kappa, A, invK = parameters[i]
        lentts = min(len(all_tts[i]), L - 1)
        popt, pcov = fit_two_scale_persistence(all_spB[i][:lentts], all_tts[i][:lentts])
        alpha, lam1, lam2 = popt
        print("fitted two scale", popt)
        alpha_err, lam1_err, lam2_err = np.sqrt(np.diag(pcov))
        all_alpha.append(alpha)
        all_alpha_err.append(alpha_err)
        all_lam1.append(lam1)
        all_lam1_err.append(lam1_err)
        all_lam2.append(lam2)
        all_lam2_err.append(lam2_err)
        all_lam0.append(-1/np.log(all_tts[i][1]))
        # axs[2].semilogy(all_spB[i][:20], all_tts[i][:20], "s", mfc="None", label=rf"{xpar[i]}")
        axs[3].semilogy(all_spB[i], all_tts[i], "s", mfc="None", label=rf"{xpar[i]}")
        axs[3].semilogy(all_spB[i], tts_two_scale(all_spB[i], alpha, lam1, lam2), "--", color=axs[2].lines[-1].get_color())

    axs[3].set_xlabel(r"$s/B$")
    axs[3].set_ylabel(r"$\left<\cos{\theta}(s)\right>$")
    axs[3].legend(title=xparam, ncol=2)

    axs[4].plot(xpar, all_lam0, "--", label=r"$\lambda_0$")
    print("all_lam0", all_lam0)
    print("all_lam1", all_lam1)
    axs[4].errorbar(xpar, all_lam1, yerr=all_lam1_err, marker="o", label=r"$\lambda_1$")
    axs[4].errorbar(xpar, all_lam2, yerr=all_lam2_err, marker="s", label=r"$\lambda_2$")
    axs[4].set_xlabel(xparam)
    axs[4].set_ylabel(r"persistance length")
    axs[4].legend()

    plt.tight_layout()
    plt.savefig(f"{folder}/obs_{xparam}{tag}.png")
    plt.close()


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


def tts_two_scale(s, alpha, lam1, lam2):
    tts = (1 - alpha) * np.exp(-s / lam1) + alpha * np.exp(-s / lam2)
    return tts


def fit_two_scale_persistence(spB, tts):
    popt, pcov = curve_fit(tts_two_scale, spB, tts)
    return popt, pcov


def plot_R_distribution(filename):
    data = np.genfromtxt(filename, delimiter=",", skip_header=4)
    X, Y, Z, R = data[:, 2], data[:, 3], data[:, 4], data[:, 7]
    Rcalc = np.sqrt(X**2 + Y**2 + Z**2)
    print("R-Rcalc", R - Rcalc)
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(211, projection="3d")
    ax2 = fig.add_subplot(212)
    # ax.plot(X/R, Y/R, Z/R, "o")
    ax.plot(X, Y, 0, ".", ms=1)
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax2.hist(X, bins=100, histtype="step", label="X")
    ax2.hist(Y, bins=100, histtype="step", label="Y")
    ax2.hist(Z, bins=100, histtype="step", label="Z")
    ax2.legend()
    # plt.savefig(filename.replace(".csv", ".png"))
    plt.show()
    plt.close()



def calc_structure_factor(X, Y, r):
    Sq = [1.0 / n for i in range(len(q))]
    for k in range(len(q)):
        Sqk = 0
        qk = q[k]
        for i in range(n - 1):
            for j in range(i + 1, n):
                Sqk += 2.0 * np.sin(qk * (i - j)) / (qk * (i - j)) / (n * n)
        Sq[k] += Sqk
    return np.array(Sq)
