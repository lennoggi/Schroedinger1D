import numpy as np
import h5py
from matplotlib import pyplot as plt
import os

# ***** Parameters *****
filename = "../Data.h5"
outdir   = "Frames"

wf_ymin = -5.5
wf_ymax =  5.5

Fwf_ymin = -0.12
Fwf_ymax =  0.12

figsize = (10., 7.)
dpi     = 200

time_x =  0.65
time_y = -0.08
# **********************


if os.path.isdir(outdir):
    raise RuntimeError(f"Directory '{outdir}' already exists. Please remove it or change the 'outdir' variable before re-running this script.")

with h5py.File(filename, "r") as f:
    t         = f["Time"][()]
    x         = f["Position"][()]
    p         = f["Momentum"][()]
    V         = f["Potential"][()]
    out_every = f["Output frequency"][()]

    NT = len(t)
    NX = len(x)

    assert len(p) == NX
    assert len(V) == NX

    assert len(out_every) == 1
    out_every = out_every[0]

    for n in range(NT):
        it = n*out_every

        # Plot the wave function and its inverse DFT (should be the same)
        re_wf   = f[f"Wave function/Real part/Iteration {it}"][()]
        re_iFwf = f[f"Inverse DFT of the wave function/Real part/Iteration {it}"][()]
        outpath = f"{outdir}/re_wf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$x$", fontsize = 12.)
        plt.ylabel(r"$\operatorname{Re}\left(\psi\right)\left(x\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(wf_ymin, wf_ymax)
        plt.plot(x,       V, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "black",     label = "Potential")
        plt.plot(x,   re_wf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
        plt.plot(x, re_iFwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()


        im_wf   = f[f"Wave function/Imaginary part/Iteration {it}"][()]
        im_iFwf = f[f"Inverse DFT of the wave function/Imaginary part/Iteration {it}"][()]
        outpath = f"{outdir}/im_wf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$x$", fontsize = 12.)
        plt.ylabel(r"$\operatorname{Im}\left(\psi\right)\left(x\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.plot(x,       V, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "black",     label = "Potential")
        plt.plot(x,   im_wf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
        plt.plot(x, im_iFwf, ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()


        outpath = f"{outdir}/mod_wf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$x$", fontsize = 12.)
        plt.ylabel(r"$\left\vert\psi\right\vert\left(x\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(0., wf_ymax)
        ##plt.ylim(wf_ymin/5., wf_ymax)
        plt.plot(x, V, ls = "-", lw = 0.5, marker = ".", markersize = 1., color = "black", label = "Potential")
        plt.plot(x, np.sqrt(  re_wf*re_wf   +   im_wf*im_wf),   ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "red",       label = "Wave function")
        plt.plot(x, np.sqrt(re_iFwf*re_iFwf + im_iFwf*im_iFwf), ls = "--", lw = 0.5, marker = ".", markersize = 0.2, color = "limegreen", label = "Inverse DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()



        # Plot the DFT of wave function
        re_Fwf   = f[f"DFT of the wave function/Real part/Iteration {it}"][()]
        outpath  = f"{outdir}/re_Fwf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$p$", fontsize = 12.)
        plt.ylabel(r"$\operatorname{Re}\left[\hat{\psi}\right]\left(p\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(Fwf_ymin, Fwf_ymax)
        plt.plot(p, re_Fwf, ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "blue",  label = "DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()


        im_Fwf   = f[f"DFT of the wave function/Imaginary part/Iteration {it}"][()]
        outpath  = f"{outdir}/im_Fwf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$p$", fontsize = 12.)
        plt.ylabel(r"$\operatorname{Im}\left[\hat{\psi}\right]\left(p\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(Fwf_ymin, Fwf_ymax)
        plt.plot(p, im_Fwf, ls = "-", lw = 0.5, marker = ".", markersize = 1., color = "blue",  label = "DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()


        outpath = f"{outdir}/mod_Fwf"
        os.makedirs(outpath, exist_ok = True)

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$p$", fontsize = 12.)
        plt.ylabel(r"$\left\vert\hat{\psi}\right\vert\left(p\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(0., Fwf_ymax)
        plt.plot(x, np.sqrt(re_Fwf*re_Fwf + im_Fwf*im_Fwf), ls = "-",  lw = 0.5, marker = ".", markersize = 1.,  color = "blue",  label = "DFT of the wave function")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{outpath}/it_{n:06d}.png")
        plt.close()

        print(f"Iteration {it}/{NT*out_every}, time {t[n]:.2e} done")
