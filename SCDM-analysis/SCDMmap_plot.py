import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm          # for colorbar zero / normalization
from time import perf_counter
from os import path
import sys

scdm_list = ["SCDMarray_BCL11A.npy","SCDMarray_bcl-cluster1-2-3-4.npy","SCDMarray_bcl-cluster1-2.npy",
             "SCDMarray_bcl-cluster1.npy","SCDMarray_bcl-cluster2.npy","SCDMarray_bcl-cluster3-4-5.npy",
             "SCDMarray_bcl-cluster3.npy","SCDMarray_bcl-cluster4.npy","SCDMarray_bcl-cluster5-6-7-8.npy",
             "SCDMarray_bcl-cluster5.npy","SCDMarray_bcl-cluster6-7-8.npy","SCDMarray_bcl-cluster6.npy",
             "SCDMarray_bcl-cluster7.npy","SCDMarray_bcl-cluster8.npy","SCDMarray_bcl-nls.npy","SCDMarray_bcl-total.npy"]

scdm_dict = {i:np.load(i) for i in scdm_list}
for outfile,scdm in scdm_dict.items():
    PSTY = "seaborn-v0_8"
    plt.style.use(PSTY)

    SMALL_SIZE = 10
    MEDIUM_SIZE = 13
    LARGE_SIZE = 14

    plt.rc('font', size=MEDIUM_SIZE)        # controls default text sizes
    plt.rc('axes', titlesize=LARGE_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

    fig, ax = plt.subplots(figsize=(7,7))
    cnorm = TwoSlopeNorm(0, vmin=-60, vmax=60)
    img = ax.imshow(scdm, cmap="seismic", interpolation=None, norm=cnorm)
    fig.colorbar(img, ax=ax, shrink=0.75)
    N=835 # bcl11a
    diag = list(range(N))
    ax.plot(diag,diag,"k-")
    ax.set_xlim(0,N-1)
    ax.set_ylim(N-1,0)
    ax.set_xlabel(r"$j$")
    ax.set_ylabel(r"$i$")
    full_seqname = outfile[:-4]
    ax.set_title(r"SCDM" + " : " + full_seqname.upper())
    plt.tight_layout()
    SaveDir = '.'
    if SaveDir:
        outfig = outfile[:-4]
        outfig += ".pdf"
        outfpath = path.join(SaveDir,outfig)
        print(f"Saving SCDM plot as:\t'{outfpath}'\n")
        fig.savefig(outfpath)
    plt.show()
    plt.close()
