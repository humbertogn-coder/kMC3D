#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sei_chemistry_figure_cfg.py
═══════════════════════════
Reproduces the user's "SEI CHEMISTRY" sparkline figure: for each of the
9 reaction-flux variables, shows the normalized trend vs k0 and the
percent change low→high k0, with Mann-Whitney significance.

Works from a per-cycle ML-style CSV that has columns: k0 and the 9 vars.
By default it reads the reaction_tracking + ce_capacity merged data for
the sweep, but you can point it at any CSV with --csv.

Run:
  python sei_chemistry_figure_cfg.py --config config_plating
  python sei_chemistry_figure_cfg.py --config config_plating --csv my.csv
"""
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib, matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import mannwhitneyu
from scipy.signal import savgol_filter
import kmc_io as io
warnings.filterwarnings("ignore")

# The 9 variables and display names (exactly the user's set)
VARS = [
    ('salt_IDP1',           'Salt decomp. (group 1)'),
    ('salt_IDP2',           'Salt decomp. (group 2)'),
    ('solv_IDP1',           'Solvent decomp. (group 1)'),
    ('solv_IDP3',           'Solvent decomp. (group 3)'),
    ('LiF Formation',       'LiF formation'),
    ('Li2O Formation',      'Li₂O formation'),
    ('RxnPlatingSEI',       'Rxn plating → SEI'),
    ('Li Stripping',        'Li stripping'),
    ('Plating on Li metal', 'Plating on Li metal'),
]
VAR_COLORS = {
    'salt_IDP1':'#2ca25f','salt_IDP2':'#1b9e77','solv_IDP1':'#66a61e',
    'solv_IDP3':'#7570b3','LiF Formation':'#e6ab02','Li2O Formation':'#e7298a',
    'RxnPlatingSEI':'#1f78b4','Li Stripping':'#d95f02',
    'Plating on Li metal':'#e31a1c',
}

# Mapping from our reaction_tracking columns to the user's var names,
# in case the sweep CSV uses different column labels.
COLUMN_ALIASES = {
    'salt_IDP1': ['salt_IDP1', 'RxnFSI'],            # FSI reduction ~ group1
    'salt_IDP2': ['salt_IDP2', 'RxnSFO'],            # salt decomp ~ group2
    'solv_IDP1': ['solv_IDP1', 'RxnSOL'],
    'solv_IDP3': ['solv_IDP3', 'RxnF5D'],
    'LiF Formation': ['LiF Formation', 'F_xyz'],
    'Li2O Formation': ['Li2O Formation', 'O_xyz'],
    'RxnPlatingSEI': ['RxnPlatingSEI'],
    'Li Stripping': ['Li Stripping', 'RxnStripping'],
    'Plating on Li metal': ['Plating on Li metal', 'RxnPlating'],
}


def resolve_columns(df):
    """Map each desired var to an available column in df."""
    resolved = {}
    for var, _ in VARS:
        for cand in COLUMN_ALIASES.get(var, [var]):
            if cand in df.columns:
                resolved[var] = cand
                break
    return resolved


def get_trend(df, col, k0_vals):
    return np.array([df[np.isclose(df['k0'], k)][col].mean() for k in k0_vals])

def smooth1d(arr, w=7):
    s = pd.Series(arr).interpolate(limit_direction='both')
    try:
        wl = min(w, len(s)//2*2-1)
        return savgol_filter(s.values, wl, 2) if wl >= 3 else s.values
    except Exception:
        return s.values

def norm01(arr):
    mn, mx = np.nanmin(arr), np.nanmax(arr)
    return (arr-mn)/(mx-mn) if mx > mn else arr*0+0.5

def lighten(hex_color, factor=0.78):
    h=hex_color.lstrip('#')
    r,g,b=int(h[0:2],16),int(h[2:4],16),int(h[4:6],16)
    r=int(r+(255-r)*factor); g=int(g+(255-g)*factor); b=int(b+(255-b)*factor)
    return f'#{r:02x}{g:02x}{b:02x}'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",required=True)
    ap.add_argument("--csv",default=None,
                    help="CSV with per-cycle k0+variables (default: build "
                         "from reaction_tracking + ce_capacity in out_dir)")
    ap.add_argument("--k0-low",type=float,default=None,
                    help="threshold: cases with k_p <= this form the LOW group")
    ap.add_argument("--k0-high",type=float,default=None,
                    help="threshold: cases with k_p >= this form the HIGH group")
    ap.add_argument("--kp-min",type=float,default=None,
                    help="only keep cases with k_p >= this (e.g. the reference) "
                         "to make a zoomed figure of a sub-range")
    ap.add_argument("--kp-max",type=float,default=None,
                    help="only keep cases with k_p <= this")
    ap.add_argument("--suffix",default="",
                    help="appended to the output filename, e.g. _ref_to_high")
    args=ap.parse_args()
    cfg=io.load_config(args.config)
    out=cfg["out_dir"]

    # Load data
    if args.csv:
        df=pd.read_csv(args.csv)
    else:
        # Merge reaction_tracking + ce_capacity on (k0, cycle)
        rt=out/"reaction_tracking_per_cycle.csv"
        cc=out/"ce_capacity_per_cycle.csv"
        if not rt.exists():
            print(f"Need {rt}; run reaction_tracking_cfg.py first."); return
        df=pd.read_csv(rt)
        if cc.exists():
            dcc=pd.read_csv(cc)
            df=df.merge(dcc[["k0","cycle","CE","capacity_C"]],
                        on=["k0","cycle"],how="left")
    print(f"Loaded {len(df)} rows, k0 values: {df['k0'].nunique()}")

    resolved=resolve_columns(df)
    missing=[v for v,_ in VARS if v not in resolved]
    if missing:
        print(f"WARNING: no column found for: {missing}")
        print("  (figure will skip those rows)")

    # Optional zoom on a k_p sub-range (used for the reference -> highest figure)
    if args.kp_min is not None:
        df = df[df['k0'] >= args.kp_min * (1 - 1e-9)]
    if args.kp_max is not None:
        df = df[df['k0'] <= args.kp_max * (1 + 1e-9)]
    if df.empty:
        print("No rows left after the k_p range filter."); return
    k0_vals=sorted(df['k0'].unique())
    if len(k0_vals)<3:
        print("Not enough k0 values for a trend figure."); return
    log_k0=np.log10(np.array(k0_vals))
    x_norm=(log_k0-log_k0.min())/(log_k0.max()-log_k0.min())

    # Low/high thresholds default: lowest and highest third of k0 range
    k0_low = args.k0_low if args.k0_low else np.percentile(k0_vals, 25)
    k0_high= args.k0_high if args.k0_high else np.percentile(k0_vals, 75)
    low_df=df[df['k0']<=k0_low]; high_df=df[df['k0']>=k0_high]
    print(f"Low k0 ≤ {k0_low:g}: {low_df['k0'].nunique()} conditions; "
          f"High k0 ≥ {k0_high:g}: {high_df['k0'].nunique()} conditions")

    plt.rcParams.update({'font.family':'serif',
        'font.serif':['DejaVu Serif'],'font.size':9.5,'axes.linewidth':0.8,
        'xtick.direction':'in','figure.dpi':150,'savefig.dpi':300})

    avail=[(v,l) for v,l in VARS if v in resolved]
    N=len(avail)
    if N==0: print("No variables available."); return
    TITLE_H=0.15; BOTTOM_H=0.065
    fig=plt.figure(figsize=(8.8,7.4),facecolor='white')
    gs=gridspec.GridSpec(N,3,left=0.01,right=0.975,
        top=1.0-TITLE_H-0.01,bottom=BOTTOM_H,hspace=0.20,wspace=0.06,
        width_ratios=[2.6,5.2,2.2])
    hdr_y=1.0-TITLE_H-0.005
    fig.text(0.385,hdr_y,r'Trend vs $k_p$  (normalized)',ha='center',va='top',
             fontsize=8.5,color='#444',fontstyle='italic')
    fig.text(0.845,hdr_y,r'$\Delta$%  low $\to$ high $k_p$',ha='center',va='top',
             fontsize=8.5,color='#444',fontstyle='italic')
    fig.add_artist(matplotlib.lines.Line2D([0.02,0.97],[hdr_y-0.020,hdr_y-0.020],
        transform=fig.transFigure,color='#bbb',lw=0.7))
    k0_ticks=[k0_vals[0],k0_vals[len(k0_vals)//4],k0_vals[len(k0_vals)//2],
              k0_vals[3*len(k0_vals)//4],k0_vals[-1]]
    xt=[(np.log10(k)-log_k0.min())/(log_k0.max()-log_k0.min()) for k in k0_ticks]

    for i,(var,label) in enumerate(avail):
        col=resolved[var]; is_last=(i==N-1)
        lc=VAR_COLORS[var]; fc=lighten(lc)
        raw=get_trend(df,col,k0_vals); sm=smooth1d(raw)
        tn=norm01(sm); rn=norm01(raw)
        lm=low_df[col].mean(); hm=high_df[col].mean()
        pct=(hm-lm)/abs(lm)*100 if abs(lm)>1e-12 else 0.0
        try:
            _,p=mannwhitneyu(low_df[col].dropna(),high_df[col].dropna(),
                             alternative='two-sided')
            sig='***' if p<0.001 else('**' if p<0.01 else('*' if p<0.05 else 'ns'))
        except Exception:
            sig='ns'
        ax_l=fig.add_subplot(gs[i,0]); ax_l.axis('off')
        ax_l.text(1.0,0.5,label,transform=ax_l.transAxes,ha='right',va='center',
                  fontsize=9.2,color='#111')
        ax_s=fig.add_subplot(gs[i,1])
        ax_s.fill_between(x_norm,0,tn,color=fc,alpha=0.45,zorder=1)
        ax_s.scatter(x_norm,rn,s=8,color=lc,alpha=0.75,zorder=3,linewidths=0)
        ax_s.plot(x_norm,tn,color=lc,lw=1.8,zorder=4)
        ax_s.axvline(x_norm[0],color='#aaa',lw=0.5,ls=':'); ax_s.axvline(x_norm[-1],color='#aaa',lw=0.5,ls=':')
        ax_s.axhline(0,color='#ddd',lw=0.4); ax_s.axhline(1,color='#ddd',lw=0.4,ls='--')
        ax_s.set_xlim(-0.02,1.02); ax_s.set_ylim(-0.18,1.28)
        for sp in ax_s.spines.values(): sp.set_visible(False)
        ax_s.set_yticks([])
        if i==0:
            ax_s.text(0.01,1.18,'Low $k_p$',transform=ax_s.transAxes,fontsize=7.5,
                      color='#666',fontstyle='italic',va='bottom')
            ax_s.text(0.99,1.18,'High $k_p$',transform=ax_s.transAxes,fontsize=7.5,
                      color='#666',fontstyle='italic',va='bottom',ha='right')
        if is_last:
            ax_s.spines['bottom'].set_visible(True); ax_s.spines['bottom'].set_linewidth(0.7)
            ax_s.set_xticks(xt); ax_s.set_xticklabels([f'{k:g}' for k in k0_ticks],fontsize=8.2)
            ax_s.tick_params(axis='x',length=3,width=0.7,pad=2)
            ax_s.set_xlabel(r'$k_p$',fontsize=10,labelpad=4)
        else: ax_s.set_xticks([])
        ax_b=fig.add_subplot(gs[i,2])
        ax_b.barh(0,pct,color=lc,alpha=0.86,height=0.50,zorder=3)
        ax_b.axvline(0,color='#444',lw=0.8,zorder=5)
        pad=abs(pct)*0.07+0.6; xpos=pct+(pad if pct>=0 else -pad)
        ax_b.text(xpos,0,f'{pct:+.1f}% {sig}',va='center',
                  ha='left' if pct>=0 else 'right',fontsize=8.2,color='#111')
        lim=max(36,abs(pct)*1.3)
        ax_b.set_xlim(-lim,lim); ax_b.set_ylim(-0.55,0.55)
        for sp in ax_b.spines.values(): sp.set_visible(False)
        ax_b.set_yticks([])
        if is_last:
            ax_b.spines['bottom'].set_visible(True); ax_b.spines['bottom'].set_linewidth(0.7)
            ax_b.set_xlabel(r'$\Delta$%',fontsize=10,labelpad=4)
            ax_b.tick_params(axis='x',length=3,width=0.7,labelsize=8,pad=2)
        else: ax_b.set_xticks([])

    tc='#238b45'; t0=1.0-0.006
    fig.text(0.50,t0,'S E I   C H E M I S T R Y',ha='center',va='top',fontsize=8.0,
             color=tc,fontweight='bold')
    fig.text(0.50,t0-0.030,rf'Effect of $k_p$ on Reaction Fluxes ({cfg["sweep_name"]})',
             ha='center',va='top',fontsize=12.5,fontweight='bold',color='#111')
    fig.text(0.50,t0-0.063,rf'$k_p$ range {k0_vals[0]:g} to {k0_vals[-1]:g}  ·  '
             rf'low $k_p \leq {k0_low:g}$, high $k_p \geq {k0_high:g}$  ·  '
             rf'$\Delta$% = (high $-$ low)/|low|',
             ha='center',va='top',fontsize=8,color='#555',fontstyle='italic')
    sep_y=t0-0.082
    fig.add_artist(matplotlib.lines.Line2D([0.02,0.97],[sep_y,sep_y],
        transform=fig.transFigure,color=tc,lw=1.3))
    png=out/f'fig_sei_chemistry_multicolor{args.suffix}.png'
    fig.savefig(png,dpi=300,bbox_inches='tight',facecolor='white'); plt.close(fig)
    print(f"Saved {png}")


if __name__ == "__main__":
    main()
