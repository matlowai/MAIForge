#!/usr/bin/env python3
"""Offline numeric change ledger, inventory mapping, reproducibility gates.
Usage: report_v3.py FIGURES_DIR EVIDENCE_DIR
Only consumes real generated JSON; never changes materials or masks failures.
"""
import hashlib,json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1];fig=Path(sys.argv[1]);ev=Path(sys.argv[2])
c=json.loads((root/'data/catalog_v3.json').read_text());old=json.loads((root/'data/catalog_v2.json').read_text())
def write(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
def flatten(v,path=''):
 out={}
 if isinstance(v,dict):
  for k,x in v.items():out.update(flatten(x,path+'/'+k))
 elif isinstance(v,list):
  for i,x in enumerate(v):
   key=x.get('id',x.get('mineral',str(i))) if isinstance(x,dict) else str(i)
   out.update(flatten(x,path+'/'+key))
 elif isinstance(v,(float,int)) and not isinstance(v,bool):out[path]=v
 return out
def changes(a,b):
 aa,bb=flatten(a),flatten(b)
 return [dict(path=k,before=aa.get(k),after=bb.get(k)) for k in sorted(aa.keys()|bb.keys()) if aa.get(k)!=bb.get(k)]
core=changes(old,c['core']);write(ev/'core_numeric_changes.json',core)
base={r['id']:r for r in c['core']['rocks']}
variants={r['id']:changes(base[r['lithofacies']],r['spec']) for r in c['realizations']}
write(ev/'realization_numeric_changes.json',variants)
resolved=json.loads((fig/'resolved_numbers_v3.json').read_text());specs={r['id']:r['spec'] for r in c['realizations']}
body_changes={r['id']:dict(minerals=changes(c['core']['minerals'],r['resolved_minerals']),spec=changes(specs[r['id']],r['resolved_spec'])) for r in resolved}
write(ev/'resolved_numeric_changes.json',body_changes)
extensions={k:v for k,v in c.items() if k not in ['core','realizations']}
extensions['realizations']=[{k:v for k,v in r.items() if k!='spec'} for r in c['realizations']]
added=changes({},extensions);write(ev/'assemblage_numeric_additions.json',added)
# Source drift check covers CURRENT producer, not just a hard-coded test list.
inv=root.parent/'pf-inventory/src/petrology.rs';source=inv.read_text();block=source.split('pub fn mineral_density(')[1].split('/// V2')[0]
ids=set(re.findall(r'"([a-z_]+)"\s*=>',block));mapped={x['inventory_id'] for x in c['inventory_mapping']}
missing=sorted(ids-mapped);assert not missing,missing
write(ev/'inventory_coverage.json',dict(source=str(inv.resolve()),sha256=hashlib.sha256(inv.read_bytes()).hexdigest(),producer_ids=sorted(ids),mapping_ids=sorted(mapped),missing=missing))
lines=['# Inventory to pf-rockspec v3 mapping','','Canonical mineral IDs and dry-solid mass fractions are pf-rockspec IDs. Snapshot source SHA256: `'+hashlib.sha256(inv.read_bytes()).hexdigest()+'` (`pf-inventory/src/petrology.rs`). The source was read only. The verification script checks its density-match vocabulary against this mapping.','','`w_i = phi_i rho_i / sum(phi_j rho_j)`. Supply and record every producer density in kg/m3. `phi` is solid volume fraction, not bulk volume including pores. A supplied species split is a dry-mass split within that group. No producer density is replaced silently.','','| Producer ID (solid volume) | Supported pf-rockspec ID(s) | Required choice |','|---|---|---|']
for m in c['inventory_mapping']:lines.append('| '+m['inventory_id']+' | '+', '.join(m['rockspec_ids'])+' | '+('Explicit body-admitted species/optical population split' if m['requires_species_resolution'] else 'Identity mapping; source density receipt')+' |')
lines+=['','Group candidates are the supported catalogue subset, not a claim of mineralogical synonymy. Clinopyroxene must be explicitly qualified as augite; amphibole as hornblende. Clay must resolve to the admitted illite/smectite/kaolinite fractions. Silicic glass must resolve its absorption population (pale felsic glass versus the existing dark glass proxy). Missing evidence is an error, not an average.','','The next arc joins the actual inventory modes and provenance to body realizations. This mapping does not overwrite those modes or change pf-inventory. pf-fracture consumes the same canonical mass modes, mechanics in Pa and Pa sqrt(m), lengths in m, and friction angles in radians.']
(root/'INVENTORY_MAPPING_V3.md').write_text('\n'.join(lines)+'\n')
rows=json.loads((fig/'samples_v3.json').read_text())['samples'];lum=lambda x:sum(a*b for a,b in zip(x,[.2126,.7152,.0722]))
optics=[]
for r in rows:
 out=r['mean_2m'];v2=r['mean_v2_reference'];optics.append(dict(id=r['id'],v3_Y=lum(out['base_reflectance_linear']),v2_family_Y=None if v2 is None else lum(v2['base_reflectance_linear']),v3_f0=out['f0_linear'],v3_density=out['dry_bulk_density_kg_m3'],v2_density=None if v2 is None else v2['dry_bulk_density_kg_m3']))
write(ev/'optical_comparison_v3.json',optics)
write(ev/'surface_numeric_changes.json',{r['id']:changes(r['mean_v2_reference'] or {},r['mean_2m']) for r in rows})
text=['# Every changed number, v3','','All numeric priors, intervals, modes, grains, porosities, joints, mechanics and derived densities are enumerated by JSON path. `null` means absent in the compared version; it is not zero. Numeric array indexes are RGB channels or interval endpoints. Categorical/semantic changes are in MODEL_V3.md.','','- `core_numeric_changes.json`: v2 catalogue -> v3 core (including new families/minerals).','- `realization_numeric_changes.json`: each selected realization spec relative to its v3 facies reference.','- `resolved_numeric_changes.json`: body parameter effects on minerals and recomputed spec density.','- `assemblage_numeric_additions.json`: all new metadata, fabric, mechanics ratios, parameters, plume and overlay priors.', '- `surface_numeric_changes.json`: every numeric surface-output difference from its v2 family.', '- `optical_comparison_v3.json`: summary of every realized 2 m dry mean against its v2 family.','- Full canonical new values: `data/catalog_v3.json`, `figures/resolved_numbers_v3.json`.','',f'Core changed/added/removed numeric fields: {len(core)}. Realization spec differences: {sum(map(len,variants.values()))}. Resolved body differences: {sum(len(x["minerals"])+len(x["spec"]) for x in body_changes.values())}. Added assemblage numeric fields: {len(added)}.','','## 2 m fresh dry family comparison','','F8 optical fixes are intentionally inherited wherever their minerals occur. Corrected family joints also change the existing optical slope/escape calculation. The sandstone/arkose mode correction substantially brightens their v3 references; no tint is fitted to preserve the old number. These are physically motivated input changes, not a calibration pass. The optical law is still v2, with its coarse-grain dark bias. v1/v2 results themselves remain byte-identical.','','| Family | v2 Y | v3 reference Y | Change |','|---|---:|---:|---:|']
for r in optics:
 if r['id'].endswith('.reference') and r['v2_family_Y'] is not None:text.append(f"| {r['id'].split('.')[0]} | {r['v2_family_Y']:.9f} | {r['v3_Y']:.9f} | {r['v3_Y']-r['v2_family_Y']:+.9f} |")
(ev.parent/'CHANGED_NUMBERS.md').write_text('\n'.join(text)+'\n')
# Every original data/source file is checked; lib.rs is the intentional version dispatch addition.
baseline=json.loads((ev/'baseline_files.json').read_text());changed=[];unchanged=[]
for file,sha in baseline.items():
 p=Path(file)
 (unchanged if hashlib.sha256(p.read_bytes()).hexdigest()==sha else changed).append(file)
assert changed==['worldfoundry/crates/pf-rockspec/src/lib.rs'],changed
before=ev/'before';after=ev/'after';matches=[]
for p in sorted(before.iterdir()):
 if p.is_file():
  q=after/p.name;assert q.exists() and p.read_bytes()==q.read_bytes(),p.name
  matches.append(dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
write(ev/'keyoff_identity.json',dict(original_unchanged_files=unchanged,intentional_source_change=changed,byte_identical_exports=matches))
print(json.dumps(dict(core_numeric_changes=len(core),realizations=len(rows),inventory_ids=len(ids),original_unchanged=len(unchanged),identical_exports=len(matches))))
