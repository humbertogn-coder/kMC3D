#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
porosity_profiles_cfg.py
════════════════════════
SEI porosity profile + SEI thickness growth per cycle, faithfully
following the user's two original scripts:
  - "Comparative porosity profiles ... Gaussian density field + z-warp"
  - "SEI Thickness vs Cycle comparison"

Key method details (copied exactly from those scripts):
  - Radii: F=1.60, O=1.81, N=1.85, S=2.22, Li=1.73, F5D=4.00, SFO=3.63
  - Scale factors computed from the ALL field (Li+SEI), applied to both
  - SEI thickness from the 2-98% CDF of the lateral solid fraction in z
  - Porosity = 1 - solid.mean(axis=(0,1))  (lateral mean per z slice)
  - Post-60 Å stabilization of the porosity tail (optional, configurable)

Outputs to CFG['out_dir']:
  porosity_profiles.csv             (z_rel, porosity_raw, porosity_corrected)
  sei_growth_per_cycle.csv          (k0, cycle, SEI_thickness_A)
  sei_morphology_ML_features.csv    one row per k0
  fig_porosity_profile_high_vs_low.png
  fig_sei_growth_high_vs_low.png
  fig_sei_thickness_vs_k0.png

Run:  python porosity_profiles_cfg.py --config config_plating
"""
import argparse, warnings, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import kmc_io as io
warnings.filterwarnings("ignore")

# ── Model constants (from user's porosity script) ───────────────────────
SEI_SPECIES = {"F","O","N","S","F5D","SFO"}
LI_SPECIES  = "Li"
IGNORE      = {"ETH","SOL","FSI"}
ANODE_Z_WIDTH = 20.2
ANODE_HALF    = ANODE_Z_WIDTH / 2.0
VOXEL, GRID_PAD = 0.8, 2.0
SIGMA_FRAC, CUTOFF_FRAC = 0.45, 2.5
DENSITY_THR = 0.12
WARP_CAP, WARP_SMOOTH, WARP_MIN_SOLID_VOX = (0.7, 2.5), 5, 10

RADIUS = {"Li":1.73,"F":1.60,"O":1.81,"N":1.85,"S":2.22,"F5D":4.00,"SFO":3.63}
DEFAULT_RADIUS = 2.0
def _r3(r): return (4.0/3.0)*math.pi*r**3
V_REAL = {"Li":_r3(1.73),"F":_r3(1.60),"O":_r3(1.81),"N":_r3(1.85),
          "S":_r3(2.22),"F5D":_r3(4.00)*1.30,"SFO":_r3(3.63)*1.25}
V_REAL_DEFAULT = _r3(2.0)

# Thickness CDF percentiles
THICKNESS_Q_LOW, THICKNESS_Q_HIGH = 0.02, 0.98
ACTIVE_SLICE_THR = 0.003

# Post-60Å stabilization (set ENABLE_TAIL_FIX=False to disable)
ENABLE_TAIL_FIX = True
TRANSITION_Z_A = 60.0
PROFILE_SMOOTH_SIGMA, TAIL_SMOOTH_SIGMA, TAIL_BLEND = 2.0, 2.0, 0.55
REF_WINDOW_LEFT_A, REF_SLOPE_CLIP = 12.0, 0.0015
FLUCT_SIGMA_MULT, MIN_FLUCT_ABS, MAX_EXTRA_UP = 1.25, 0.020, 0.055
LOWER_BAND_FACTOR, RESIDUAL_KEEP = 0.35, 0.35


def _make_kernel(radius, voxel):
    sigma = max(1e-6, SIGMA_FRAC*float(radius)); r_cut = CUTOFF_FRAC*sigma
    half = int(np.ceil(r_cut/voxel)); ax = np.arange(-half, half+1)*voxel
    X,Y,Z = np.meshgrid(ax,ax,ax,indexing="ij")
    R = np.sqrt(X**2+Y**2+Z**2)
    K = np.exp(-(R**2)/(2.0*sigma**2)); K[R>r_cut]=0.0
    return K.astype(np.float32)


def _grid(coords, voxel, pad):
    mn = coords.min(axis=0)-pad; mx = coords.max(axis=0)+pad
    dims = np.maximum(np.ceil((mx-mn)/voxel).astype(int)+1, 3)
    return mn, dims


def _stamp(field, K, idx):
    ix,iy,iz=idx; kx,ky,kz=K.shape; hx,hy,hz=kx//2,ky//2,kz//2
    x0,x1=ix-hx,ix+hx+1; y0,y1=iy-hy,iy+hy+1; z0,z1=iz-hz,iz+hz+1
    fx0,fx1=max(0,x0),min(field.shape[0],x1)
    fy0,fy1=max(0,y0),min(field.shape[1],y1)
    fz0,fz1=max(0,z0),min(field.shape[2],z1)
    if fx0>=fx1 or fy0>=fy1 or fz0>=fz1: return
    field[fx0:fx1,fy0:fy1,fz0:fz1]+=K[fx0-x0:kx-(x1-fx1),
                                      fy0-y0:ky-(y1-fy1),fz0-z0:kz-(z1-fz1)]


def build_density_fields(coords, species, voxel):
    mn, dims = _grid(coords, voxel, GRID_PAD)
    fields = {"SEI":np.zeros(dims,np.float32),"LI":np.zeros(dims,np.float32)}
    cache = {}
    for p, sp in zip(coords, species):
        r = RADIUS.get(sp, DEFAULT_RADIUS)
        if r not in cache: cache[r] = _make_kernel(r, voxel)
        idx = np.floor((p-mn)/voxel).astype(int)
        if np.any(idx<0) or np.any(idx>=dims): continue
        if sp == LI_SPECIES: _stamp(fields["LI"], cache[r], tuple(idx))
        elif sp in SEI_SPECIES: _stamp(fields["SEI"], cache[r], tuple(idx))
    fields["ALL"] = fields["SEI"] + fields["LI"]
    meta = {"origin":mn,"dims":tuple(int(v) for v in dims),"voxel":float(voxel)}
    return fields, meta


def moving_average(x, w):
    if w is None or w<=1: return x
    w=int(w);
    if w%2==0: w+=1
    pad=w//2
    return np.convolve(np.pad(x,(pad,pad),mode="edge"),np.ones(w)/w,mode="valid")


def compute_slice_scale_factors(species, coords, meta, fields):
    voxel=float(meta["voxel"]); oz=float(meta["origin"][2])
    nz=fields["ALL"].shape[2]
    solid=(fields["ALL"]>=DENSITY_THR).sum(axis=(0,1)).astype(float)
    z=coords[:,2]; iz=np.floor((z-oz)/voxel).astype(int)
    m=(iz>=0)&(iz<nz); iz=iz[m]; sp=np.asarray(species)[m]
    Vt=np.zeros(nz)
    for s,k in zip(sp,iz):
        if s in IGNORE: continue
        Vt[k]+=V_REAL.get(s,V_REAL_DEFAULT)
    Vf=solid*(voxel**3)
    scale=np.ones(nz)
    for k in range(nz):
        scale[k]=1.0 if solid[k]<WARP_MIN_SOLID_VOX else Vt[k]/max(1e-12,Vf[k])
    scale=np.clip(scale,WARP_CAP[0],WARP_CAP[1])
    return moving_average(scale,WARP_SMOOTH)


def warp_field_z(field, sf):
    nz=field.shape[2]; out=[]
    for k in range(nz):
        sk=float(sf[k]); nrep=int(np.floor(sk)); frac=sk-nrep
        sl=field[:,:,k]
        for _ in range(max(1,nrep)): out.append(sl)
        if frac>1e-6: out.append((1-frac)*sl+frac*field[:,:,min(nz-1,k+1)])
    return np.stack(out,axis=2).astype(field.dtype)


def warp_fields(fields, meta, sf):
    f_sei=warp_field_z(fields["SEI"],sf); f_li=warp_field_z(fields["LI"],sf)
    meta2=dict(meta); meta2["dims"]=(f_sei.shape[0],f_sei.shape[1],f_sei.shape[2])
    return {"SEI":f_sei,"LI":f_li,"ALL":f_sei+f_li}, meta2


def sei_thickness_from_field(field_sei, meta, up_abs_z):
    """Thickness from 2-98% CDF of lateral solid fraction in z."""
    voxel=float(meta["voxel"]); oz=float(meta["origin"][2])
    sol=(field_sei>=DENSITY_THR)
    z_profile=sol.mean(axis=(0,1))
    z_abs=oz+np.arange(z_profile.size)*voxel
    w=np.clip(z_profile,0.0,None)
    if w.sum()<1e-12: return 0.0
    cdf=np.cumsum(w)/np.sum(w)
    z_high=float(np.interp(THICKNESS_Q_HIGH,cdf,z_abs))
    return max(0.0, z_high-float(up_abs_z))


def stabilize_tail(z_rel, por_raw):
    """Stabilize porosity after TRANSITION_Z_A (user's correction)."""
    out=por_raw.copy()
    if len(out)<8 or not ENABLE_TAIL_FIX: return out
    idx=np.searchsorted(z_rel,TRANSITION_Z_A)
    if idx>=len(out)-2: return out
    smooth=gaussian_filter1d(por_raw,sigma=PROFILE_SMOOTH_SIGMA)
    z0=TRANSITION_Z_A
    win=(z_rel>=max(0.0,z0-REF_WINDOW_LEFT_A))&(z_rel<=z0)
    if win.sum()<4:
        i0=max(0,idx-6); i1=max(i0+3,idx+1); z_win=z_rel[i0:i1]; y_win=smooth[i0:i1]
    else:
        z_win=z_rel[win]; y_win=smooth[win]
    baseline=float(np.median(y_win))
    slope=np.polyfit(z_win,y_win,1)[0] if (len(z_win)>=2 and z_win.max()-z_win.min()>1e-8) else 0.0
    slope=float(np.clip(slope,-REF_SLOPE_CLIP,REF_SLOPE_CLIP))
    z_tail=z_rel[idx:]; ref=baseline+slope*(z_tail-z0)
    allow=min(max(MIN_FLUCT_ABS,FLUCT_SIGMA_MULT*float(np.std(y_win))),MAX_EXTRA_UP)
    upper=ref+allow; lower=ref-LOWER_BAND_FACTOR*allow
    tail_raw=por_raw[idx:]
    tail_smooth=gaussian_filter1d(tail_raw,sigma=TAIL_SMOOTH_SIGMA)
    tail_mix=(1-TAIL_BLEND)*tail_raw+TAIL_BLEND*tail_smooth
    high=tail_mix>upper
    tail_mix[high]=upper[high]+RESIDUAL_KEEP*(tail_mix[high]-upper[high])
    out[idx:]=np.clip(tail_mix,lower,upper)
    return out


def extract_porosity_profile(field, meta, up_abs_z):
    voxel=float(meta["voxel"]); oz=float(meta["origin"][2]); nz=field.shape[2]
    solid=(field>=DENSITY_THR).mean(axis=(0,1))
    por_raw=1.0-solid
    z_abs=oz+np.arange(nz)*voxel
    iz_anode=int(np.clip(np.searchsorted(z_abs,float(up_abs_z)),0,nz-1))
    # top via active slices
    active=solid>ACTIVE_SLICE_THR; idx_a=np.where(active)[0]
    if idx_a.size==0: return None
    w=np.clip(solid,0,None)
    cdf=np.cumsum(w)/max(np.sum(w),1e-12)
    z_top=float(np.interp(THICKNESS_Q_HIGH,cdf,z_abs))
    iz_top=int(np.clip(np.searchsorted(z_abs,z_top),0,nz-1))
    if iz_top<=iz_anode: return None
    z_rel=z_abs[iz_anode:iz_top+1]-float(up_abs_z)
    por_trim=por_raw[iz_anode:iz_top+1]
    por_corr=stabilize_tail(z_rel,por_trim)
    return {"z_rel":z_rel,"por_raw":por_trim,"por_corr":por_corr}


def process_case(case_key, k0, cfg):
    base=cfg["base_dir"]/case_key
    traj=base/cfg["traj_subpath"]; log=io.resolve_log(base, cfg)
    if not traj.is_dir() or log is None:
        print(f"  [{case_key}] missing; skip", flush=True); return None, None
    calib=traj/"kmc-coords-0.xyz"
    if not calib.is_file():
        print(f"  [{case_key}] no calib; skip", flush=True); return None, None
    fr0=io.parse_xyz_last_frame(calib)
    co0=io.minimal_image(fr0["coords"],fr0["cell"]); sp0=fr0["species"]
    z_li=co0[sp0==LI_SPECIES,2]
    up_abs_z=(float(np.median(z_li))+ANODE_HALF) if z_li.size else 0.0

    df_log=io.parse_data2excel(log)
    s2c=io.build_step_to_cycle(df_log)
    cyc2xyz=io.map_cycles_to_xyz(io.list_xyz_files(traj), s2c)
    cycles=list(range(cfg["cycle_step_reactions"],cfg["max_cycle"]+1,
                       cfg["cycle_step_reactions"]))
    g_rows, p_rows = [], []
    db_cyc = cfg.get("li_buried_cycle_for_db",110)
    for cyc in cycles:
        if cyc not in cyc2xyz: continue
        fr=io.parse_xyz_last_frame(cyc2xyz[cyc][1])
        co=io.minimal_image(fr["coords"],fr["cell"]); sp=fr["species"]
        keep=(~np.isin(sp,list(IGNORE)))&(co[:,2]>=up_abs_z)
        if keep.sum()<5: continue
        coords,species=co[keep],sp[keep]
        fields,meta=build_density_fields(coords,species,VOXEL)
        sf=compute_slice_scale_factors(species,coords,meta,fields)
        fields_w,meta_w=warp_fields(fields,meta,sf)
        thick=sei_thickness_from_field(fields_w["SEI"],meta_w,up_abs_z)
        g_rows.append({"case_key":case_key,"k0":k0,"cycle":cyc,
                       "SEI_thickness_A":thick})
        if cyc==cycles[-1] or cyc==db_cyc:
            prof=extract_porosity_profile(fields_w["SEI"],meta_w,up_abs_z)
            if prof is not None:
                for zr,pr,pc in zip(prof["z_rel"],prof["por_raw"],prof["por_corr"]):
                    p_rows.append({"case_key":case_key,"k0":k0,"cycle":cyc,
                                   "z_rel_A":float(zr),"porosity_raw":float(pr),
                                   "porosity_corrected":float(pc)})
    if not g_rows: return None, None
    dg=pd.DataFrame(g_rows); dp=pd.DataFrame(p_rows)
    print(f"  [{case_key}] final SEI thickness = {dg['SEI_thickness_A'].iloc[-1]:.1f} A",
          flush=True)
    return dg, dp


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True)
    args=ap.parse_args(); cfg=io.load_config(args.config)
    out=cfg["out_dir"]; out.mkdir(parents=True,exist_ok=True)
    print(f"Porosity & growth — {cfg['sweep_name']}")
    g_parts,p_parts=[],[]
    for case_key,k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──",flush=True)
        dg,dp=process_case(case_key,k0,cfg)
        if dg is not None: g_parts.append(dg)
        if dp is not None and not dp.empty: p_parts.append(dp)
    if not g_parts: print("No data."); return
    df_growth=pd.concat(g_parts,ignore_index=True)
    df_growth.to_csv(out/"sei_growth_per_cycle.csv",index=False)
    df_profile=pd.concat(p_parts,ignore_index=True) if p_parts else pd.DataFrame()
    if not df_profile.empty:
        df_profile.to_csv(out/"porosity_profiles.csv",index=False)
    ml=(df_growth.groupby(["case_key","k0"])
        .agg(SEI_thickness_A_mean=("SEI_thickness_A","mean"),
             SEI_thickness_A_final=("SEI_thickness_A","last"))
        .reset_index())
    ml.to_csv(out/"sei_morphology_ML_features.csv",index=False)
    print(f"\nSaved {len(df_growth)} growth rows, {len(ml)} k0 rows")

    plt.rcParams.update(cfg["plot_style"])
    CASES = cfg["comparison_cases"]   # (key, k0, color, label)

    # Porosity profile high vs low (corrected)
    if not df_profile.empty:
        fig,ax=plt.subplots(figsize=(8.0,5.0))
        for key,k0,color,label in CASES:
            g=io.select_by_k0(df_profile,k0)
            if g.empty: continue
            cyc=g["cycle"].max(); gc=g[g["cycle"]==cyc].sort_values("z_rel_A")
            ax.plot(gc["z_rel_A"],gc["porosity_raw"],ls="--",lw=1.0,
                    alpha=0.30,color=color)
            ax.plot(gc["z_rel_A"],gc["porosity_corrected"],
                    lw=2.7 if label=="reference" else 2.1,
                    color=color,label=f"{cfg['rate_symbol_plain']} = {k0:g} ({label})")
        if ENABLE_TAIL_FIX:
            ax.axvline(TRANSITION_Z_A,color="gray",ls=":",lw=1.2,alpha=0.8)
        ax.set_xlabel("Distance from anode surface (Å)")
        ax.set_ylabel("SEI porosity"); ax.set_title("Comparative SEI porosity profiles")
        ax.set_ylim(-0.02,1.02); ax.set_xlim(left=0.0)
        ax.grid(alpha=0.25); ax.legend(frameon=True)
        plt.tight_layout()
        fig.savefig(out/"fig_porosity_profile_high_vs_low.png",bbox_inches="tight")
        plt.close()

    # SEI thickness vs cycle
    fig,ax=plt.subplots(figsize=(8.0,5.2))
    for key,k0,color,label in CASES:
        g=io.select_by_k0(df_growth,k0)
        if g.empty: continue
        ax.plot(g["cycle"],g["SEI_thickness_A"],
                lw=2.6 if label=="reference" else 2.0,color=color,
                label=f"{cfg['rate_symbol_plain']} = {k0:g} ({label})")
    ax.set_xlabel("Cycle"); ax.set_ylabel("SEI Thickness (Å)")
    ax.set_title("Comparative SEI Thickness vs Cycle")
    ax.grid(alpha=0.25); ax.legend(frameon=True)
    plt.tight_layout()
    fig.savefig(out/"fig_sei_growth_high_vs_low.png",bbox_inches="tight"); plt.close()

    # SEI thickness vs k0
    fig,ax=plt.subplots(figsize=(8.0,5.2))
    fin=(df_growth.sort_values("cycle").groupby("k0").last().reset_index().sort_values("k0"))
    ax.plot(fin["k0"],fin["SEI_thickness_A"],"o-",lw=2.0,color="#34699A",markersize=7)
    ax.set_xscale("log"); ax.set_xlabel(f"{cfg['rate_symbol_plain']} (plating rate)")
    ax.set_ylabel("SEI thickness, final cycle (Å)")
    ax.set_title(f"SEI thickness vs {cfg['rate_symbol_plain']}"); ax.grid(alpha=0.25,which="both",lw=0.5)
    plt.tight_layout()
    fig.savefig(out/"fig_sei_thickness_vs_k0.png",bbox_inches="tight"); plt.close()
    print("Saved figures.")


if __name__ == "__main__":
    main()
