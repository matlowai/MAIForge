#!/usr/bin/env python3
"""Build source-labelled all-key sanity table from exported evaluator outputs."""
import json,argparse
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('evidence',type=Path);a=ap.parse_args();p=a.evidence
v1=json.loads((p/'swatches_v1.json').read_text())['samples'];v2=json.loads((p/'swatches_v2.json').read_text())['samples']
def Y(x):return sum(a*b for a,b in zip(x['output']['base_reflectance_linear'],[.2126,.7152,.0722]))
def find(rows,key,f):return next(x for x in rows if x['input']['material_ref']['spec_id']==key and x['input']['footprint_m']==f and x['input']['surface']['surface_id'].endswith('fresh_dry'))
bands={k:(float(lo),float(hi)) for k,lo,hi in (x.split() for x in '''shale .08 .20
mudstone .10 .22
siltstone .15 .30
greywacke .10 .20
sandstone_quartz .25 .50
arkose .25 .40
conglomerate .20 .40
limestone .30 .60
dolostone .30 .50
chalk .60 .85
chert .10 .35
coal .03 .06
rock_salt .50 .80
slate .06 .15
phyllite .10 .20
schist .10 .25
gneiss .15 .35
quartzite .30 .60
marble .50 .90
rhyolite .20 .40
pumice .40 .60
sand .30 .50
silt_clay .15 .35
till .15 .30
regolith .15 .35
basalt_tholeiitic .05 .10
granite .25 .35
obsidian .008 .03'''.splitlines())}
lit={k:tuple(float(v)/100 for v in [b,g,r]) for k,b,g,r in (x.split() for x in '''shale 17.90 17.87 20.45
mudstone 29.47 38.54 47.45
limestone 47.68 48.59 47.87
sandstone_quartz 35.51 38.78 35.80
conglomerate 24.21 24.25 24.67
granite 30.54 27.43 27.03
basalt_tholeiitic 16.17 15.46 16.27
quartzite 38.27 45.08 46.74
marble 78.46 89.48 89.40
gneiss 19.24 17.78 18.22'''.splitlines())}
# Basalt alkali is left unknown: no false type-specific source mapping.
keys=list(dict.fromkeys(x['input']['material_ref']['spec_id'] for x in v2))
rows=[];table=['# All-rock v1/v2 sanity comparison','','Metric: linear pre-light diffuse Y. V2 primary comparison is a dry fresh 2 m outcrop footprint. V1 had identical values at all three footprints. Also show v2 0.1 m and 50 m. No displayed-sRGB or total hemispherical albedo equivalence is asserted.','','Literature: [Lee, Jeong & Kim 2021, Table 3](https://www.engeojournal.org/articles/xml/b2nd/). The cited interval is the min/max of the three published visible spectral-band averages for a specimen, NOT a lithology population or uncertainty interval. Four sampled locations; no population uncertainty supplied. Ten matched/proxy specimens; 33 keys have UNKNOWN numeric source bands. Basalt and sandstone are class proxies, not specific subtype measurements. A source-span miss is not automatically a physics failure.','','| Rock | v1 Y | v2 0.1 m | v2 2 m | v2 50 m | Cited spectral span | Source comparison | Integrator band | Brief comparison |','|---|---:|---:|---:|---:|---|---|---|---|']
for k in keys:
 old=Y(find(v1,k,.1));values=[Y(find(v2,k,f)) for f in [.1,2.,50.]];value=values[1];band=bands.get(k);source=lit.get(k);sb=[min(source),max(source)] if source else None
 status='unknown' if band is None else 'pass' if band[0]<=value<=band[1] else 'miss'
 lstatus='unknown' if sb is None else 'inside' if sb[0]<=value<=sb[1] else 'outside'
 row=dict(rock=k,v1_Y=old,v2_Y_01m=values[0],v2_Y_2m=value,v2_Y_50m=values[2],integrator_band=band,integrator_comparison=status,cited_band=sb,cited_band_kind='spectral_band_span_not_uncertainty' if sb else 'unknown',literature_BGR=source,literature_comparison=lstatus,reference='Lee et al. 2021 Table 3; DOI 10.9720/kseg.2021.4.479' if source else None,source_proxy=k in ['sandstone_quartz','basalt_tholeiitic'])
 rows.append(row);fmt=lambda x:'unknown' if x is None else f'{x[0]:.3f}-{x[1]:.3f}'
 table.append(f'| {k} | {old:.3f} | {values[0]:.3f} | {value:.3f} | {values[2]:.3f} | {fmt(sb)} | {lstatus} | {fmt(band)} | {status} |')
table+=['','Measured-source disagreement with brief: the selected mudstone is .385 visible mean versus .10-.22 brief; basalt .160 versus .05-.10. Preserve those disagreements. Marble has .858 mean; granite .283. These do not justify replacing one conditional composition with another solely to hit a band.','','Obsidian .008-.03 is an arc-1 preservation envelope chosen in this arc. The earlier brief supplied no numeric obsidian band. Marble has an additional protected lower bound .60 tested at all three footprints. Other missing brief bands are unknown.']
(p/'sanity_table_v2.md').write_text('\n'.join(table)+'\n');(p/'sanity_v2.json').write_text(json.dumps(dict(samples=rows),indent=2)+'\n')
print(json.dumps(dict(rows=len(rows),brief_pass=sum(r['integrator_comparison']=='pass' for r in rows),brief_miss=[r['rock'] for r in rows if r['integrator_comparison']=='miss'],brief_unknown=sum(r['integrator_comparison']=='unknown' for r in rows),literature_known=sum(r['cited_band'] is not None for r in rows))))
