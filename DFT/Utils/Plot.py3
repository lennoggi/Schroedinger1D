import numpy as np
import h5py
from matplotlib import pyplot as plt

# ***** Parameters *****
filename = "../Data.h5"
figname  = "DFT.pdf"
# **********************

with h5py.File(filename, "r") as f:
    x           = f["Position"][()]
    p           = f["Momentum"][()]
    re_wf       = f["Wave function/Real part"][()]
    im_wf       = f["Wave function/Imaginary part"][()]
    re_Fwf      = f["DFT of the wave function/Real part"][()]
    im_Fwf      = f["DFT of the wave function/Imaginary part"][()]
    re_true_Fwf = f["Analytical DFT of the wave function/Real part"][()]
    im_true_Fwf = f["Analytical DFT of the wave function/Imaginary part"][()]
    re_iFwf     = f["Inverse DFT of the wave function/Real part"][()]
    im_iFwf     = f["Inverse DFT of the wave function/Imaginary part"][()]

    N = len(x)
    assert len(      re_wf == N)
    assert len(      im_wf == N)
    assert len(     re_Fwf == N)
    assert len(     im_Fwf == N)
    assert len(re_true_Fwf == N)
    assert len(im_true_Fwf == N)
    assert len(    re_iFwf == N)
    assert len(    im_iFwf == N)


    # Plot the wave function and its inverse DFT (should be the same)
    plt.figure()
    plt.xlabel("$x$", fontsize = 12.)
    plt.ylabel(r"$\operatorname{Re}\left(\psi\right)\left(x\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(x,   re_wf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
    plt.plot(x, re_iFwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
    plt.legend()
    plt.tight_layout()
    plt.savefig("re_wf.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("$x$", fontsize = 12.)
    plt.ylabel(r"$\operatorname{Im}\left(\psi\right)\left(x\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(x,   im_wf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
    plt.plot(x, im_iFwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
    plt.legend()
    plt.tight_layout()
    plt.savefig("im_wf.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("$x$", fontsize = 12.)
    plt.ylabel(r"$\left\vert\psi\right\vert\left(x\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(x, np.sqrt(  re_wf*re_wf   +   im_wf*im_wf),   ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
    plt.plot(x, np.sqrt(re_iFwf*re_iFwf + im_iFwf*im_iFwf), ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
    plt.legend()
    plt.tight_layout()
    plt.savefig("mod_wf.pdf")
    plt.close()


    # Plot the DFT of the wave function
    plt.figure()
    plt.xlabel("$p$", fontsize = 12.)
    plt.ylabel(r"$\operatorname{Re}\left[\hat{\psi}\right]\left(p\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(p,      re_Fwf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "blue", label = "Numerical DFT")
    plt.plot(p, re_true_Fwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "gold", label = "Analytical DFT")
    plt.legend()
    plt.tight_layout()
    plt.savefig("re_Fwf.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("$p$", fontsize = 12.)
    plt.ylabel(r"$\operatorname{Im}\left[\hat{\psi}\right]\left(p\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(p,      im_Fwf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "blue", label = "Numerical DFT")
    plt.plot(p, im_true_Fwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "gold", label = "Analytical DFT")
    plt.legend()
    plt.tight_layout()
    plt.savefig("im_Fwf.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("$p$", fontsize = 12.)
    plt.ylabel(r"$\left\vert\hat{\psi}\right\vert\left(p\right)$", fontsize = 12.)
    plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
    plt.plot(x, np.sqrt(re_Fwf*re_Fwf + im_Fwf*im_Fwf),
             ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "blue", label = "Numerical DFT")
    plt.plot(x, np.sqrt(re_true_Fwf*re_true_Fwf + im_true_Fwf*im_true_Fwf),
             ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "gold", label = "Analytical DFT")
    plt.legend()
    plt.tight_layout()
    plt.savefig("mod_Fwf.pdf")
    plt.close()
