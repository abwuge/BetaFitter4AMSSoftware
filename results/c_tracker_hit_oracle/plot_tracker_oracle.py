#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE=Path(__file__).resolve().parent
OUT=HERE/"figures"
OUT.mkdir(exist_ok=True)
P6=["#5790fc","#f89c20","#e42536","#964a8b","#9c9ca1","#7a21dd"]

def setup():
    plt.rcParams.update({
        "font.size":11,"axes.labelsize":12,"axes.titlesize":13,
        "legend.fontsize":9.5,"xtick.labelsize":10.5,"ytick.labelsize":10.5,
        "figure.figsize":(7.4,4.9),"figure.dpi":150,"savefig.dpi":220,
        "axes.linewidth":0.9,"lines.linewidth":2.0,
        "pdf.fonttype":42,
    })

def clean(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out",length=4,width=0.9)
    ax.grid(axis="y",alpha=0.18,linewidth=0.7)

def save(fig,name,pdf):
    fig.tight_layout()
    fig.savefig(str(OUT/(name+".png")),bbox_inches="tight")
    fig.savefig(str(OUT/(name+".pdf")),bbox_inches="tight")
    pdf.savefig(fig,bbox_inches="tight")
    plt.close(fig)

def main():
    setup()
    strict=json.load(open(HERE/"strict_results.json"))
    transfer=json.load(open(HERE/"strict_transfer_results.json"))
    unsafe=json.load(open(HERE/"transfer_results.json"))
    with PdfPages(str(OUT/"tracker_oracle_summary.pdf")) as pdf:
        stages=["strict_multifit","strict_fitdiag","strict_geometry","strict_morph"]
        labels=["multi-fit"," + residual/fit diag"," + hit geometry"," + COG morphology"]
        fig,ax=plt.subplots()
        x=np.arange(len(stages))
        y=np.array([strict["bootstrap"][s]["overall"]["median"] for s in stages])
        lo=np.array([strict["bootstrap"][s]["overall"]["low95"] for s in stages])
        hi=np.array([strict["bootstrap"][s]["overall"]["high95"] for s in stages])
        for i,s in enumerate(stages):
            ax.errorbar(x[i],y[i],yerr=[[y[i]-lo[i]],[hi[i]-y[i]]],fmt="o",capsize=4,
                        color=P6[i],markersize=7,label=labels[i])
        ax.set_xticks(x); ax.set_xticklabels(labels,rotation=12,ha="right")
        ax.set_ylabel(r"$W_{68}(\Delta\ln R)$ reduction vs raw [%]")
        ax.set_title("Strict tracker-only MC oracle: recoverable rigidity resolution")
        ax.set_ylim(0,max(10,1.15*np.max(hi)))
        clean(ax); save(fig,"01_strict_overall_gain",pdf)

        fig,ax=plt.subplots()
        centers=np.array([2.75,3.5,4.5,5.5])
        for i,s in enumerate(stages):
            bins=strict["bootstrap"][s]["bins"]
            yy=np.array([b["median"] for b in bins])
            ll=np.array([b["low95"] for b in bins])
            hh=np.array([b["high95"] for b in bins])
            ax.errorbar(centers,yy,yerr=[yy-ll,hh-yy],marker="o",capsize=3,
                        color=P6[i],label=labels[i])
        ax.set_xticks(centers); ax.set_xticklabels(["2.5–3","3–4","4–5","5–6"])
        ax.set_xlabel("reconstructed rigidity [GV]")
        ax.set_ylabel(r"$W_{68}(\Delta\ln R)$ reduction vs raw [%]")
        ax.set_title("Recoverable rigidity information by rigidity bin")
        ax.legend(frameon=False,ncol=2)
        clean(ax); save(fig,"02_strict_gain_vs_rigidity",pdf)

        fig,ax=plt.subplots()
        tests=["C12","C13"]; xx=np.arange(2)
        for i,stage in enumerate(["strictgeom_c12train","strictgeom_c13train"]):
            b=transfer["bootstrap"][stage]
            yy=np.array([b[t]["median"] for t in tests])
            ll=np.array([b[t]["low95"] for t in tests])
            hh=np.array([b[t]["high95"] for t in tests])
            ax.errorbar(xx+(i-0.5)*0.10,yy,yerr=[yy-ll,hh-yy],fmt="o",capsize=4,
                        markersize=7,color=P6[i],label=("trained on C12" if i==0 else "trained on C13"))
        ax.set_xticks(xx); ax.set_xticklabels([r"test $^{12}$C",r"test $^{13}$C"])
        ax.set_ylabel(r"$W_{68}(\Delta\ln R)$ reduction vs raw [%]")
        ax.set_title("Strict hit-geometry correction transfers across isotopes")
        ax.legend(frameon=False)
        clean(ax); save(fig,"03_strict_cross_isotope_transfer",pdf)

        fig,ax=plt.subplots()
        items=[
            ("strictgeom_c12train",transfer["results"],"strict: train C12",P6[0]),
            ("strictgeom_c13train",transfer["results"],"strict: train C13",P6[1]),
            ("charge_c12train",unsafe["results"],"charge: train C12",P6[2]),
            ("charge_c13train",unsafe["results"],"charge: train C13",P6[3]),
        ]
        offsets=[-0.18,-0.06,0.06,0.18]
        for (stage,src,label,color),off in zip(items,offsets):
            vals=[src[stage]["by_isotope"][t]["corrected_median"] if "by_isotope" in src[stage] else src[stage][t]["median"] for t in tests]
            ax.plot(xx+off,vals,"o",color=color,markersize=7,label=label)
        ax.axhline(0,color=P6[4],linestyle=":",linewidth=1.3)
        ax.set_xticks(xx); ax.set_xticklabels([r"test $^{12}$C",r"test $^{13}$C"])
        ax.set_ylabel(r"median corrected $\Delta\ln R$")
        ax.set_title("Charge shortcut narrows widths by creating isotope-dependent scale shifts")
        ax.legend(frameon=False,ncol=2)
        clean(ax); save(fig,"04_charge_shortcut_scale_bias",pdf)
    print("wrote",OUT/"tracker_oracle_summary.pdf")

if __name__=="__main__":
    main()
