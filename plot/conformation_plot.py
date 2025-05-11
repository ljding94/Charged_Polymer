import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib import rc
import os


def get_obs_data(folder, L, kappa, A, invK):
    finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
    filename = f"{folder}/obs_{finfo}.csv"
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)

    return data


def get_obs_data_per_invK(folder, L, kappa, A, invKs):
    all_R2 = []
    all_R2_err = []
    all_Rg2 = []
    all_Rg2_err = []
    all_Sq = []
    all_q = []
    all_tts = []
    all_spB = []

    for i in range(len(invKs)):
        data = get_obs_data(folder, L, kappa, A, invKs[i])
        E, E_B, E_HSY, R2, Rg2 = data[0, 5], data[0, 6], data[0, 7], data[0, 8], data[0, 9]
        E_err, E_B_err, E_HSY_err, R2_err, Rg2_err = data[1, 5], data[1, 6], data[1, 7], data[1, 8], data[1, 9]

        Sq = data[0, 10:]
        q = data[2, 10:]
        tts = data[3, 10:]
        spB = data[5, 10:]

        all_R2.append(R2 / L**2)
        all_R2_err.append(R2_err / L**2)
        all_Rg2.append(Rg2 / L**2)
        all_Rg2_err.append(Rg2_err / L**2)

        all_Sq.append(Sq)
        all_q.append(q)
        all_tts.append(tts)
        all_spB.append(spB)

    return all_R2, all_R2_err, all_Rg2, all_Rg2_err, all_Sq, all_q, all_tts, all_spB


def ax_plot_config_data(ax, folder, L, kappa, A, invK, xshift=0, yshift=0):
    # get config coordination
    finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
    filename = f"{folder}/config_{finfo}.csv"
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]

    ax.plot(x + xshift, y + yshift, "-", color="gray")


def ax_plot_multi_config_data(ax, folder, L, kappa, A, invK, xshift=0, yshift=0):
    finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
    for i in range(30):
        filename = f"{folder}/{finfo}/config_{i}.csv"
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        x, y, z = data[:, 0], data[:, 1], data[:, 2]
        X, Y, Z = x[-1], y[-1], z[-1]
        R = np.sqrt(X**2 + Y**2 + Z**2)
        # color = (0.4*X/R+0.6, (0.4*Y/R+0.6)*1.0, 0.4*Z/R+0.6)
        color = (np.abs(Y) / R, np.abs(Z) / R, np.abs(X) / R)
        ax.plot(x + xshift, y + yshift, "-", color=color, lw=0.75, solid_capstyle="round")  # , rasterized=True)
    ax.set_aspect("equal")


def plot_conformation_kappa_A(tex_lw=240.71031, ppi=72):
    print("\nplotting R2 Rg2 figure\n")
    # plot Rg and R2 versus invK(lamD) for various kappa
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax00 = plt.subplot2grid((2, 2), (0, 0))
    ax01 = plt.subplot2grid((2, 2), (0, 1), sharex=ax00, sharey=ax00)
    ax10 = plt.subplot2grid((2, 2), (1, 0), sharex=ax00)
    ax11 = plt.subplot2grid((2, 2), (1, 1), sharex=ax00, sharey=ax10)

    folder = "../data/20250115"
    L = 500
    A = 5.0
    invKs = np.arange(1.0, 10.1, 0.5)
    for kappa in [10, 30, 50]:
        all_R2, all_R2_err, all_Rg2, all_Rg2_err, all_Sq, all_q, all_tts, all_spB = get_obs_data_per_invK(folder, L, kappa, A, invKs)

        # ax00.errorbar(invKs, all_R2, yerr=np.array(all_R2_err), ls="none", marker="o", ms=3, mfc="none", lw=1, label=f"{kappa:.0f}")
        ax00.plot(invKs, all_R2, "o", ms=3, mfc="none", lw=1, label=f"{kappa:.0f}")

        # ax10.errorbar(invKs, all_Rg2, yerr=np.array(all_Rg2_err), ls="none", marker="s", ms=3, mfc="none", lw=1, label=f"{kappa:.0f}")
        ax10.plot(invKs, all_Rg2, "s", ms=3, mfc="none", lw=1, label=f"{kappa:.0f}")
    kappa = 30
    for A in [3, 5, 7]:
        all_R2, all_R2_err, all_Rg2, all_Rg2_err, all_Sq, all_q, all_tts, all_spB = get_obs_data_per_invK(folder, L, kappa, A, invKs)

        # ax01.errorbar(invKs, all_R2, yerr=all_R2_err, ls="none", marker="o", ms=3, mfc="none", lw=1, label=f"{A:.0f}")
        ax01.plot(invKs, all_R2, "o", ms=3, mfc="none", lw=1, label=f"{A:.0f}")

        # ax11.errorbar(invKs, all_Rg2, yerr=all_Rg2_err, ls="none", marker="s", ms=3, mfc="none", lw=1, label=f"{A:.0f}")
        ax11.plot(invKs, all_Rg2, "s", ms=3, mfc="none", lw=1, label=f"{A:.0f}")

    ax00.set_ylabel(r"$R_2/L^2$", fontsize=9, labelpad=0)
    # ax1.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelsize=7)
    ax00.legend(title=r"$\kappa$",loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax00.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax00.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax00.yaxis.set_major_locator(plt.MultipleLocator(0.3))
    ax00.yaxis.set_minor_locator(plt.MultipleLocator(0.15))

    ax10.set_ylabel(r"$R_g^2/L^2$", fontsize=9, labelpad=0)
    ax10.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelsize=7, labelbottom=True, labelleft=True)
    ax10.legend(title=r"$\kappa$", loc="upper left",ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax10.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax10.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax10.yaxis.set_major_locator(plt.MultipleLocator(0.02))
    ax10.yaxis.set_minor_locator(plt.MultipleLocator(0.01))

    # ax01.set_ylabel(r"$R_2/L^2$", fontsize=9, labelpad=0)
    # ax1.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=False, labelsize=7)
    ax01.legend(title=r"$A$",loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # ax11.set_ylabel(r"$R_g^2/L^2$", fontsize=9, labelpad=0)
    ax11.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelsize=7, labelbottom=True, labelleft=False)
    ax11.legend(title=r"$A$",loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # plot sample configurations
    # L, A, kappa = 500, 5, 30
    # invKs = [2, 4, 6, 8]
    # for i in range(len(invKs)):
    #    ax_plot_config_data(ax3, folder, L, kappa, A, invKs[i], yshift=i * 100)
    # ax3.set_axis_off()
    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]
    for ax in [ax00, ax01, ax10, ax11]:
        ax.text(0.8, 0.15, annotation.pop(0), transform=ax.transAxes, fontsize=9)

    plt.tight_layout(pad=0.5)
    plt.savefig("./figures/R2_Rg2.pdf", format="pdf")
    plt.savefig("./figures/R2_Rg2.png", dpi=300)
    plt.show()
    plt.close()


# TODO: get back to this, may not needed, since not much information provided
def plot_configuration(tex_lw=240.71031, ppi=72):
    print("\nplotting configuration figure\n")
    # plot Rg and R2 versus invK(lamD) for various kappa
    fig = plt.figure(figsize=(tex_lw / ppi * 2, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax0 = plt.subplot2grid((1, 2), (0, 0))
    ax0f = plt.subplot2grid((1, 2), (0, 0), zorder=-1)
    ax1 = plt.subplot2grid((1, 2), (0, 1))
    ax1f = plt.subplot2grid((1, 2), (0, 1), zorder=-1)
    shift = 700

    # kappa vs invK
    kappas = [10, 30, 50]
    invKs = [3, 5, 9]

    ax0f.set_ylim(0, len(invKs))
    ax0f.set_yticks(np.linspace(1, len(invKs), len(invKs)) - 0.5)
    ax0f.set_yticklabels([rf"{invKs[i]:.0f}" for i in range(len(invKs))])
    ax0f.set_ylabel(r"$\lambda_D$", fontsize=9, labelpad=0)

    ax0f.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax0f.set_xlim(0, len(kappas))
    ax0f.set_xticks(np.linspace(1, len(kappas), len(kappas)) - 0.5)
    ax0f.set_xticklabels([rf"{kappas[i]:.0f}" for i in range(len(kappas))])
    ax0f.set_xlabel(r"$\kappa$", fontsize=9, labelpad=0)
    ax0f.set_title(r"$A=9$", fontsize=9)
    ax0f.set_aspect("equal")

    for j in range(len(kappas)):
        for i in range(len(invKs)):
            pass
            ax_plot_multi_config_data(ax0, "../data/20250118_cfg", 500, kappas[i], 9, invKs[j], xshift=i * shift, yshift=j * shift)
    ax0.set_axis_off()
    # ax0.set_xlim(-shift*0.25, (len(invKs)-0.5) * shift)
    # A vs invK
    As = [1, 3, 7]


    ax1f.set_ylim(0, len(invKs))
    ax1f.set_yticks(np.linspace(1, len(invKs), len(invKs)) - 0.5)
    ax1f.set_yticklabels([rf"{invKs[i]:.0f}" for i in range(len(invKs))])
    ax1f.set_ylabel(r"$\lambda_D$", fontsize=9, labelpad=0)

    ax1f.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax1f.set_xlim(0, len(As))
    ax1f.set_xticks(np.linspace(1, len(As), len(As)) - 0.5)
    ax1f.set_xticklabels([rf"{As[i]:.0f}" for i in range(len(As))])
    ax1f.set_xlabel(r"$A$", fontsize=9, labelpad=0)
    ax1f.set_title(r"$\kappa=30$", fontsize=9)
    ax1f.set_aspect("equal")

    for j in range(len(As)):
        for i in range(len(invKs)):
            ax_plot_multi_config_data(ax1, "../data/20250118_cfg", 500, 30, As[i], invKs[j], xshift=i * shift, yshift=j * shift)
    ax1.set_axis_off()

    plt.tight_layout(pad=1)
    plt.savefig("./figures/configuration.pdf", format="pdf")
    plt.savefig("./figures/configuration.png", dpi=300)
    plt.show()
    plt.close()


def get_obs_tts(folder, L, kappa, A, invK):
    finfo = f"L{L}_kappa{kappa:.1f}_A{A:.1f}_invK{invK:.1f}"
    filename = f"{folder}/obs_{finfo}.csv"
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    tts = data[3, 10:]
    spB = data[5, 10:]
    return spB, tts


def tts_one_scale(s, lam0):
    tts = np.exp(-s / lam0)
    return tts


def fit_one_scale_persistence(spB, tts):
    popt, pcov = curve_fit(tts_one_scale, spB, tts)
    return popt, pcov


def tts_two_scale(s, alpha, lam1, lam2):
    tts = (1 - alpha) * np.exp(-s / lam1) + alpha * np.exp(-s / lam2)
    # tts = 1 - (1 - alpha) * s / lam1 + alpha * (np.exp(-s / lam2) - 1)
    return tts


def fit_two_scale_persistence(spB, tts):
    bounds = (0, [1.0, np.inf, np.inf])
    popt, pcov = curve_fit(tts_two_scale, spB, tts, bounds=bounds)
    return popt, pcov


def get_all_peristence_lengths(spB, tts):
    lam0, alpha, lam1, lam2 = np.nan, np.nan, np.nan, np.nan
    # one scale fit
    popt0, pcov0 = fit_one_scale_persistence(spB, tts)
    lam0 = popt0[0]
    popt1, pcov1 = fit_two_scale_persistence(spB, tts)
    alpha, lam1, lam2 = popt1[0], popt1[1], popt1[2]
    if lam1 < lam2:
        print(f"lam1 < lam2: lam1={lam1}, lam2={lam2}")
        alpha, lam1, lam2 = 1 - alpha, lam2, lam1
    return lam0, alpha, lam1, lam2


def plot_two_length_scale(tex_lw=240.71031, ppi=72):
    print("\nplotting two length scale figure\n")
    # plot Rg and R2 versus invK(lamD) for various kappa
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    # for tts, two ways of fitting
    ax00 = plt.subplot2grid((2, 2), (0, 0))
    ax01 = plt.subplot2grid((2, 2), (0, 1), sharex=ax00, sharey=ax00)

    # got lam, various kappa and A
    ax10 = plt.subplot2grid((2, 2), (1, 0))
    ax11 = plt.subplot2grid((2, 2), (1, 1), sharex=ax10, sharey=ax10)

    folder = "../data/20250115"
    L = 500
    kappa = 30
    A = 5

    # tts plot
    invKs = [1, 3, 5, 7, 9]
    for i in range(len(invKs)):
        spB, tts = get_obs_tts(folder, 500, 30, 5, invKs[i])
        ax00.plot(spB[::5], tts[::5], "o", ms=2.5, mfc="none", lw=1, label=f"{invKs[i]:.0f}")
        popt, pcov = fit_one_scale_persistence(spB, tts)
        lam0 = popt[0]
        ax00.plot(spB, tts_one_scale(spB, lam0), "-", lw=1, color=ax00.get_lines()[-1].get_color())

        ax01.plot(spB[::5], tts[::5], "o", ms=2.5, mfc="none", lw=1, label=f"{invKs[i]:.0f}")
        popt1, pcov1 = fit_two_scale_persistence(spB, tts)
        alpha, lam1, lam2 = popt1[0], popt1[1], popt1[2]
        if lam1 < lam2:
            print(f"lam1 < lam2: lam1={lam1}, lam2={lam2}")
            alpha, lam1, lam2 = 1 - alpha, lam2, lam1
        lame = lam2 / alpha
        ax01.plot(spB, tts_two_scale(spB, alpha, lam1, lam2), "-", lw=1, color=ax00.get_lines()[-1].get_color())

    ax00.set_ylabel(r"$\left<cos\theta(s) \right>$", fontsize=9, labelpad=0)
    ax00.set_xlabel(r"$s$", fontsize=9, labelpad=0)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax00.legend(title=r"$\lambda_D$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9, labelspacing=0.2)
    ax00.xaxis.set_major_locator(plt.MultipleLocator(30))
    ax00.xaxis.set_minor_locator(plt.MultipleLocator(15))
    ax00.yaxis.set_major_locator(plt.MultipleLocator(0.3))
    ax00.yaxis.set_minor_locator(plt.MultipleLocator(0.15))
    ax00.set_xlim(0, 125)

    # ax01.set_ylabel(r"$\left<cos\theta \right>(s)$", fontsize=9, labelpad=0)
    ax01.set_xlabel(r"$s$", fontsize=9, labelpad=0)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=7)
    ax01.legend(title=r"$\lambda_D$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9, labelspacing=0.2)

    # TODO: finish this
    # 1. variation of kappa
    kappas = [10, 30, 50]
    colors = ["royalblue", "tomato", "limegreen"]
    cuts = [5, 6, 9]
    invKs = np.arange(1.0, 10.1, 0.5)

    for i in range(len(kappas)):
        lam0s, alphas, lam1s, lam2s, lames = [], [], [], [], []
        for j in range(len(invKs)):
            spB, tts = get_obs_tts(folder, 500, kappas[i], 5, invKs[j])
            lam0, alpha, lam1, lam2 = get_all_peristence_lengths(spB, tts)
            lam0s.append(lam0)
            lam1s.append(lam1)
            alphas.append(alpha)
            lam2s.append(lam2)
            lames.append(lam2 / alpha)
        ax10.plot(invKs[: cuts[i] + 1], lam0s[: cuts[i] + 1], "-", ms=3, mfc="none", lw=1, color=colors[i], label=f"{kappas[i]:.0f}")
        ax10.plot(invKs[cuts[i] :], lam1s[cuts[i] :], "--", ms=3, mfc="none", lw=1, color=colors[i])
        ax10.plot(invKs[cuts[i] :], lames[cuts[i] :], ":", ms=3, mfc="none", lw=1, color=colors[i])

    ax10.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax10.set_ylabel(r"$\lambda$", fontsize=9, labelpad=0)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax10.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax10.xaxis.set_major_locator(plt.MultipleLocator(3))
    ax10.xaxis.set_minor_locator(plt.MultipleLocator(1.5))
    ax10.yaxis.set_major_locator(plt.MultipleLocator(80))
    ax10.yaxis.set_minor_locator(plt.MultipleLocator(40))
    ax10.set_ylim(-30, 290)

    # 2. variation of A
    As = [1, 5, 9]
    cuts = [9, 7, 6]
    for i in range(len(As)):
        lam0s, alphas, lam1s, lam2s, lames = [], [], [], [], []
        for j in range(len(invKs)):
            spB, tts = get_obs_tts(folder, 500, 30, As[i], invKs[j])
            lam0, alpha, lam1, lam2 = get_all_peristence_lengths(spB, tts)
            lam0s.append(lam0)
            lam1s.append(lam1)
            alphas.append(alpha)
            lam2s.append(lam2)
            lames.append(lam2 / alpha)
        ax11.plot(invKs[: cuts[i] + 1], lam0s[: cuts[i] + 1], "-", ms=3, mfc="none", lw=1, color=colors[i], label=f"{As[i]:.0f}")
        ax11.plot(invKs[cuts[i] :], lam1s[cuts[i] :], "--", ms=3, mfc="none", lw=1, color=colors[i])
        ax11.plot(invKs[cuts[i] :], lames[cuts[i] :], ":", ms=3, mfc="none", lw=1, color=colors[i])
    ax11.set_xlabel(r"$\lambda_D$", fontsize=9, labelpad=0)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=7)
    ax11.legend(title=r"$A$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]
    for ax in [ax00, ax01, ax10, ax11]:
        ax.text(0.825, 0.075, annotation.pop(0), transform=ax.transAxes, fontsize=9)

    plt.tight_layout(pad=0.1)
    plt.savefig("./figures/two_length_scale.pdf", format="pdf")
    plt.savefig("./figures/two_length_scale.png", dpi=300)
    plt.show()
    plt.close()
