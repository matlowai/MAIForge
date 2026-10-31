#!/usr/bin/env python3
"""Q24 exact receiver edges, no smoothing. Retained Q20 crop is evaluation only.
Adapted from Q20 scripts/plot.py; original hash is recorded in lane BASE.json.
"""
import argparse, json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);a=p.parse_args();run=a.run
m=json.loads((run/'metrics.json').read_text());d=json.loads((run/'domain.json').read_text());pin=json.loads((run/'inputs.json').read_text())['c7']
q20=pathlib.Path(d['q20']);sel=d['figure_selection'];extent=sel['crop_extent_deg'];center=sel['rotation_center_lon'];ny,nx=pin['shape']
wrap=lambda x:(x-center+180)%360-180+center
z=np.fromfile(pathlib.Path(pin['snapshot'])/'elevation_m.f32',dtype='<f4').reshape(ny,nx)
mask=np.zeros(nx*ny);mask[d['cells']]=1;mask[d['core']]=2;mask=mask.reshape(ny,nx)
lon=np.arange(nx)*360/nx-180+180/nx;lat=np.arange(ny)*180/ny-90+90/ny
xx=wrap(lon);ix=np.argsort(xx);xx=xx[ix];z=z[:,ix];mask=mask[:,ix]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
fig,axs=plt.subplots(1,3,figsize=(19,7),layout='constrained')
titles=['Native C7 surface pass (398–400 My)','Rejected Q20 per-cell routing','Q24 joint r8 catchment (400 My)']
paths=[q20/'coarse.csv',q20/'refined.csv',run/'joint.csv']
for ax,title,path in zip(axs,titles,paths):
    ax.imshow(z,origin='lower',extent=[xx[0]-180/nx,xx[-1]+180/nx,-90,90],cmap='gist_earth',vmin=-100,vmax=4500,alpha=.5,interpolation='nearest',aspect='auto')
    ax.contour(xx,lat,mask,levels=[.5,1.5],colors=['#888888','#c14228'],linewidths=[.8,1.2])
    data=np.loadtxt(path,delimiter=',',skiprows=1,ndmin=2);seg=data[:,1:5].reshape(-1,2,2);seg[:,:,0]=wrap(seg[:,:,0]);q=data[:,5]
    visible=(np.abs(seg[:,1,0]-seg[:,0,0])<180)&(seg[:,:,0].max(axis=1)>=extent[0])&(seg[:,:,0].min(axis=1)<=extent[1])&(seg[:,:,1].max(axis=1)>=extent[2])&(seg[:,:,1].min(axis=1)<=extent[3]);seg=seg[visible];q=q[visible]
    lc=LineCollection(seg,linewidths=np.clip(.25+.22*np.log10(q/10),.25,1.5),colors='#124e78',alpha=.95);ax.add_collection(lc)
    ax.set(xlim=extent[:2],ylim=extent[2:],xlabel='longitude (degrees)',ylabel='latitude (degrees)',title=title)
    ax.grid(alpha=.15)
    # Native cell lattice is intentionally visible as a reference, not a routing constraint.
    ax.set_xticks(np.arange(np.ceil(extent[0]*nx/360)*360/nx,extent[1],360/nx),minor=True)
    ax.set_yticks(np.arange(np.ceil(extent[2]*ny/180)*180/ny,extent[3],180/ny),minor=True)
    ax.grid(which='minor',alpha=.08)
au=m['audit'];fig.suptitle(f"Same retained Q20 crop · {m['coarse_core_cells']:,} upstream/core cells + {m['halo_cells']:,} halo · no smoothed paths\nJoint: {au['uphill_boundary_heads']} uphill boundary heads; {m['coarse_receiver_divergences']:,} receiver divergences; water residual {au['relative_residual']:.2e}",fontsize=15)
fig.supxlabel(f"Blue: Q ≥ 10 m³/s. Red outline: complete upstream catchment. Grey: halo. Native/Q20 also show surrounding catchments; Q24 shows only this domain's tagged runoff.\nRouting heads ≠ beds: {au['filled_nodes']:,} filled nodes (max {au['max_fill_m']:.1f} m), {au['uphill_boundary_bed_links']:,} uphill boundary BED links retained and disclosed.\nQ24 uses endpoint terrain + preserved native lake beds; Q20 used a pre-pass reconstruction and portal incision. Server/browser UNINTEGRATED. Nested boundary pass: {m['nested_gate_pass']}.",fontsize=10)
fig.savefig(run/'comparison.png',dpi=160);plt.close(fig)
# Detailed requested-cell view: parent and child; plots retained even if boundary gate fails.
c=m['request_cell'];cx=wrap((c%nx+.5)*360/nx-180);cy=(c//nx+.5)*180/ny-90
fig,axs=plt.subplots(1,2,figsize=(12,6),layout='constrained')
for ax,name,path in zip(axs,['Parent r8','Requested cell r32'],[run/'joint.csv',run/'nested.csv']):
    data=np.loadtxt(path,delimiter=',',skiprows=1,ndmin=2);data=data[data[:,0]==c];seg=data[:,1:5].reshape(-1,2,2);seg[:,:,0]=wrap(seg[:,:,0]);q=data[:,5]
    ax.add_collection(LineCollection(seg,linewidths=np.clip(.2+.3*np.log10(q/10),.2,1.8),colors='#124e78'))
    ax.set(xlim=[cx-180/nx,cx+180/nx],ylim=[cy-90/ny,cy+90/ny],xlabel='longitude',ylabel='latitude',title=name);ax.grid(alpha=.2)
na=m['coupled_nested'];fig.suptitle(f"Coupled nested boundary consistency: {na['boundary_consistency_pass']}\nParent head change {na['max_parent_head_change_m']:.3g} m; uphill interface heads {na['uphill_interface_heads']}; cell water residual {na['cell_relative_residual']:.2e}")
fig.supxlabel(f"Exact receiver edges; {na['changed_parent_receivers']} parent receivers corrected; {na['replacement_cycle_nodes']} cycle nodes. Full coupled product required. Server UNINTEGRATED.");fig.savefig(run/'nested.png',dpi=160);plt.close(fig)
print(run/'comparison.png');print(run/'nested.png')
