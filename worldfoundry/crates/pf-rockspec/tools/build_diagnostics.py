#!/usr/bin/env python3
"""Frozen manufactured diagnostic state, not natural-world samples."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'data/catalog_v1.json').read_text())
def reference(r): return dict(spec_id=r['id'],spec_version=c['schema_version'],realization_id='diagnostic-nominal-v1-'+r['id'],body_id='diagnostic-fixture-body')
def bulk(r):
 return dict(mineral_mode=[dict(mineral=m['mineral'],mass_fraction=m['mass_fraction']['value']) for m in r['mode']],intrinsic_porosity=r['intrinsic_porosity']['value'],interparticle_porosity=r['interparticle_porosity']['value'],saturation=0.0,temperature_k=293.15,damage=0.0,fabric_to_material=[1.,0.,0.,0.],alteration_record='manufactured nominal project-prior assemblage; no site admission')
def surface(r,history):
 s=dict(surface_id=r['id']+'/'+history,creation_event='diagnostic-cut-v1',exposure_elapsed_s=0.,forcing_record='manufactured dry fracture diagnostic',material_coordinate_m=[.123,.456,.789],normal_material=[0.,0.,1.],mechanism='fresh_fracture',wet_coverage=0.,water_film_m=0.,abrasion=0.,rms_slope_history_rad=0.,coating=None)
 if history=='fresh_wet':
  s.update(wet_coverage=1.,water_film_m=.0001,forcing_record='manufactured fully wet fracture; 100 micrometre water film')
 if history=='weathered_dry':
  s.update(creation_event='diagnostic-inherited-face',exposure_elapsed_s=None,forcing_record='explicit diagnostic coating and roughness history; no rate/time simulation',mechanism='inherited',rms_slope_history_rad=.12,abrasion=.1)
  products=r['weathering_products']; mineral=next((m for m in products if m in ('goethite','kaolinite','smectite','illite','serpentine')),None)
  if r['weathering_family']=='carbonate_solution':
   # Pure marble gets dissolution roughening, not invented rust.
   mineral='kaolinite' if r['id']!='marble' else None
  if mineral:
   s['coating']=dict(source_event='diagnostic-alteration-residue-v1',mineral_mode=[dict(mineral=mineral,mass_fraction=1.)],grain_diameter_m=.00001,thickness_m=.0001,coverage=.18 if r['weathering_family']=='carbonate_solution' else .65,extinction_m_inv=c['model']['rind_extinction_m_inv']['value'])
 return s
def case(r,history,footprint,label):
 return dict(id=label,material_ref=reference(r),grain_diameter_m=r['optical_domain_diameter_m']['value'],bulk=bulk(r),surface=surface(r,history),footprint_m=footprint)
cases=[]
for i,r in enumerate(c['rocks']):
 for j,h in enumerate(('fresh_dry','fresh_wet')):
  cases.append(case(r,h,[.1,2.,50.][(i+j)%3],f'parity-{len(cases):03d}'))
# Weathered state sampling evenly includes all families and all 3 footprint scales.
for i in list(range(29))+[33,34,35,38,42]:
 cases.append(case(c['rocks'][i],'weathered_dry',[.1,2.,50.][i%3],f'parity-{len(cases):03d}'))
for i,(rock,variant) in enumerate([('schist','cleavage_parallel'),('schist','cleavage_cross'),('granite','abraded'),('basalt_tholeiitic','altered_fresh'),('obsidian','micro_footprint'),('granite','zero_layer'),('marble','partial_wet'),('schist','rotated_fabric')]):
 r=next(r for r in c['rocks'] if r['id']==rock); x=case(r,'fresh_dry',.0001,f'parity-{len(cases):03d}');s=x['surface'];s['surface_id']+='-'+variant
 if variant.startswith('cleavage'): s['mechanism']='cleavage'; s['normal_material']=[0.,0.,1.] if variant.endswith('parallel') else [1.,0.,0.]
 elif variant=='abraded': s['abrasion']=.9;s['mechanism']='abraded'
 elif variant=='altered_fresh':
  x['bulk']['mineral_mode']=[dict(mineral='chlorite',mass_fraction=.6),dict(mineral='plagioclase',mass_fraction=.4)];x['bulk']['alteration_record']='manufactured buried hydrothermal chlorite replacement; retained on fresh cut'
 elif variant=='micro_footprint': x['footprint_m']=.000001
 elif variant=='zero_layer': s['coating']=dict(source_event='diagnostic-zero-layer',mineral_mode=[dict(mineral='goethite',mass_fraction=1.)],grain_diameter_m=.00001,thickness_m=0.,coverage=1.,extinction_m_inv=10000.)
 elif variant=='partial_wet': s['wet_coverage']=.5;s['water_film_m']=.00005
 elif variant=='rotated_fabric': x['bulk']['fabric_to_material']=[.7071067811865476,0.,.7071067811865476,0.];s['mechanism']='cleavage'
 cases.append(x)
assert len(cases)==128
(ROOT/'data/diagnostic_inputs_v1.json').write_text(json.dumps(dict(schema='pf-rockspec-diagnostics/1',qualification='manufactured contract fixtures, not measured or seeded world samples',samples=cases),indent=2)+'\n')
swatches=[case(r,h,f,f'{r["id"]}/{h}/{f}') for r in c['rocks'] for f in [.1,2.,50.] for h in ['fresh_dry','fresh_wet','weathered_dry']]
(ROOT/'data/swatch_inputs_v1.json').write_text(json.dumps(dict(schema='pf-rockspec-diagnostics/1',qualification='manufactured conditional prior swatches',samples=swatches),indent=2)+'\n')
print(f'{len(cases)} frozen diagnostics; {len(swatches)} contact sheet cells')
