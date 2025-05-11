#!/opt/homebrew/bin/python3
from conformation_plot import *
from ML_Sq_plot import *


def main():

    # 1. mechanical response: R2, Rg, and configurations
    #plot_conformation_kappa_A()

    # maybe not really need this
    plot_configuration()

    # 2. two length scale analysis: two way of fitting tangent correlation and variation versus invK(lam_D) for kappa and L, maybe add A later?
    #plot_two_length_scale()

    # 3. variation of Sq, independently tuning invK, kappa and A
    #plot_Sq_example()

    # 4.1
    #plot_SVD_basis()

    # 4. plot SVD analysis on data feature distribution

    #plot_SVD_feature()

    # 5. plot learnability analysis (nearest-neighbor distance) by slice
    #plot_SVD_NND()

    # plot LML
    #plot_LML()

    # 6. plot GPR learning results example
    #plot_GPR()


if __name__ == "__main__":
    main()
