#!/usr/bin/env python3
"""Offline v3 conditional catalogue. Never writes v1/v2. No empirical calibration.
All new numerical priors cite RS3: remembered mineralogy/mechanics, broad support.
"""
import copy, json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
old=json.loads((ROOT/'data/catalog_v2.json').read_text()); c=copy.deepcopy(old)
def q(v,lo,hi,unit='1',cond='conditional ambient dry specimen; offline remembered prior, not measurement'):
 return dict(value=v,interval=[lo,hi],unit=unit,interval_kind='plausible_prior_support',distribution='unspecified_joint_prior',reference='RS3',evidence_kind='project_prior',support_scale='body conditional realization; not site-admitted',conditions=cond,covariance_group='v3_closed_assemblage')
def est(v,unit='1',cond='conditional family prior; offline memory, not calibration'):
 return q(v,v*.5,v*1.5,unit,cond)
def mode(d):
 assert abs(sum(d.values())-1)<1e-9
 return [dict(mineral=k,mass_fraction=q(v,max(0,v*.4),min(1,v*1.6),cond='dry solid mass fraction, joint simplex constraint')) for k,v in d.items()]
c['schema_version']='pf-rockspec/3'; c['numerical_model']='assemblage+connected-km/3'
c['references']['RS3']=dict(citation='CX-MAT-ASSEMBLAGE 2026-10-25; design review material review F2/F4/F5/F7/F8; MODEL_V3.md',qualification='Offline recalled mineralogy and engineering ranges; not extracted measurements; no spectral calibration or site admission',verification='offline_memory_unverified_project_prior')
md={m['id']:m for m in c['minerals']}
# F8: corrected aggregate priors, not a change of the F3 mixing law.
md['magnetite']['extinction_k']=[q(.8,.3,1.3,cond='effective visible extinction; n=2.4 gives F0=0.213, remembered polished reflectance ~0.20') for _ in range(3)]
for k,rgb in [('plagioclase',[.74,.75,.76]),('olivine',[.36,.43,.20])]:
 md[k]['diffuse_reflectance']=[q(v,v*.55,min(.95,v*1.4),'linear_reflectance','clean aggregate at 1 mm, fresh grains; unverified optical prior') for v in rgb]
# id,rho,Mohs,n,kRGB,cleavage,angles,quality,habit,lustre,RGB,role
rows=[
 ('galena',7500,2.5,2.2,[3.,3.,3.],3,[90],'perfect','subconchoidal','metallic',[.035]*3,'ore'),
 ('sphalerite',4000,3.75,2.37,None,6,[60,120],'perfect','conchoidal','vitreous',[.25,.13,.035],'ore'),
 ('cassiterite',6950,6.5,2.0,None,2,[90],'imperfect','uneven','vitreous',[.10,.06,.03],'ore'),
 ('chromite',4700,5.5,2.2,[.7]*3,0,[],'none','uneven','sub_metallic',[.025]*3,'ore'),
 ('native_gold',19300,2.75,1.1,[2.8,2.2,1.6],0,[],'none','hackly','metallic',[.025,.018,.006],'ore'),
 ('uraninite',9800,5.5,2.1,[.9]*3,0,[],'none','conchoidal','sub_metallic',[.02]*3,'ore'),
 ('fluorite',3180,4,1.434,None,4,[70.5,109.5],'perfect','subconchoidal','vitreous',[.62,.58,.72],'gangue'),
 ('barite',4500,3.25,1.64,None,3,[90,78,102],'good','uneven','vitreous',[.76,.74,.69],'gangue'),
 ('graphite',2260,1.5,2.6,[1.4]*3,1,[],'perfect','flaky','metallic',[.025]*3,'ore'),
 ('gibbsite',2420,2.75,1.58,None,1,[],'perfect','uneven','pearly',[.75,.72,.65],'ore'),
 ('boehmite',3030,3.5,1.65,None,1,[],'perfect','uneven','earthy',[.69,.66,.59],'ore'),
 ('diaspore',3400,6.75,1.72,None,1,[],'perfect','uneven','vitreous',[.65,.62,.54],'ore'),
 ('bornite',5100,3,2.5,[2.6,2.3,2.0],0,[],'indistinct','uneven','metallic',[.045,.018,.025],'ore'),
 ('chalcocite',5600,2.75,2.4,[2.2]*3,0,[],'indistinct','conchoidal','metallic',[.025]*3,'ore'),
 ('native_silver',10500,2.75,1.2,[3.5]*3,0,[],'none','hackly','metallic',[.025]*3,'ore'),
 ('native_copper',8960,2.75,1.2,[3.1,2.4,2.1],0,[],'none','hackly','metallic',[.03,.015,.01],'ore'),
 ('molybdenite',4730,1.25,2.5,[2.5]*3,1,[],'perfect','flaky','metallic',[.035]*3,'ore'),
 ('siderite',3960,4,1.75,None,3,[73,107],'perfect','uneven','vitreous',[.40,.31,.19],'ore'),
 ('malachite',4000,3.75,1.80,None,1,[],'perfect','uneven','vitreous',[.04,.34,.12],'secondary'),
 ('azurite',3770,3.75,1.80,None,1,[],'perfect','conchoidal','vitreous',[.03,.10,.36],'secondary'),
 ('scheelite',6060,4.75,1.92,None,2,[90],'distinct','uneven','vitreous',[.72,.68,.53],'ore'),
 ('beryl',2750,7.75,1.58,None,1,[],'poor','conchoidal','vitreous',[.49,.64,.56],'accessory'),
 ('spodumene',3180,6.75,1.67,None,2,[87,93],'perfect','uneven','vitreous',[.58,.66,.51],'ore'),
]
for key,rho,h,n,k,d,angles,quality,habit,lustre,rgb,role in rows:
 m=copy.deepcopy(md['quartz']); m.update(id=key,density_kg_m3=q(rho,rho*.94,rho*1.08,'kg/m3'),mohs=q(h,max(1,h-.5),min(10,h+.5),'Mohs'),refractive_index=q(n,max(1,n-.15),n+.2),extinction_k=None if k is None else [est(x,cond='effective RGB extinction surrogate; no measured dispersion') for x in k],cleavage=dict(directions=q(d,d,d,'count'),angles_deg=[q(a,a-3,a+3,'deg') for a in angles],quality=quality,reference='RS3'),fracture_habit=habit,lustre=lustre,diffuse_reflectance=[q(v,max(.001,v*.4),min(.95,v*1.5),'linear_reflectance') for v in rgb],alteration_products=[],plagioclase_an_fraction=None)
 c['minerals'].append(m);md[key]=m
roles={r[0]:r[-1] for r in rows}
for k in ['magnetite','hematite','goethite','ilmenite','pyrite','chalcopyrite','pyrolusite']:roles[k]='ore'
for k in ['apatite','garnet','cordierite']:roles[k]='accessory'
for k in ['kaolinite','smectite','illite','chlorite','epidote','serpentine','talc']:roles[k]='alteration'
meta={k:dict(mineral_class=roles.get(k,'rock_forming'),liberation_grain_size_m=q(.0003,.00001,.01,'m','unresolved intergrowth liberation prior; body texture must constrain it; not a recovery fraction')) for k in md}
rd={r['id']:r for r in c['rocks']}
new=[('pegmatite','granite',dict(quartz=.27,k_feldspar=.40,plagioclase=.20,muscovite=.10,beryl=.02,spodumene=.01),.03,.01,130),('skarn','hornfels',dict(garnet=.35,augite=.30,epidote=.15,calcite=.15,magnetite=.05),.002,.02,160),('banded_iron_formation','chert',dict(quartz=.50,hematite=.30,magnetite=.20),.0001,.02,180),('laterite','mudstone',dict(goethite=.35,hematite=.15,gibbsite=.25,boehmite=.10,kaolinite=.15),.0002,.25,12),('mylonite','gneiss',dict(quartz=.35,k_feldspar=.25,plagioclase=.25,muscovite=.10,biotite=.05),.00005,.01,140),('gossan','mudstone',dict(goethite=.55,hematite=.15,quartz=.25,kaolinite=.05),.0003,.30,15)]
for key,base,m,g,p,u in new:
 r=copy.deepcopy(rd[base]);r['id']=key;r['mode']=mode(m);r['intrinsic_porosity']=q(p,p*.5,min(.8,p*2));r['vesicularity']=q(0,0,0);r['optical_domain_diameter_m']=est(g,'m')
 for name,v in [('d10_m',g*.2),('d50_m',g),('d90_m',g*4)]:r['grain'][name]=est(v,'m')
 r['mechanical']['ucs_pa']=dict(status='known_prior',estimate=est(u*1e6,'Pa'));r['weathering_products']=list(m);r['weathering_family']='oxidation_residual' if key in ['laterite','gossan'] else 'silicate_hydrolysis'
 r['fabric']='foliated' if key=='mylonite' else 'bedded' if key=='banded_iron_formation' else 'massive'
 r['structure_v2']['kind']='banded' if key in ['mylonite','banded_iron_formation'] else 'crystalline';r['structure_v2']['wavelength_m']=est(.02 if key=='banded_iron_formation' else g,'m');r['structure_v2']['optical_boundary_efficiency']=q(1,.2,1)
 c['rocks'].append(r);rd[key]=r
# F7: ratio chain + independent toughness by family, scalar reference direction.
# tensile/UCS, E/UCS (both Pa), KIC MPa sqrt(m), friction deg
chains={'crystalline_silicate':(.06,350,1.8,38),'carbonate':(.07,650,1.2,33),'cemented_clastic':(.07,400,.8,32),'argillaceous':(.04,250,.5,21),'evaporite':(.08,1000,.35,27),'glass':(.05,580,.8,30),'porous_volcanic':(.08,300,.4,32),'foliated':(.05,400,1.2,28),'organic':(.04,200,.25,25),'residual':(.06,250,.3,28),'unconsolidated':(0,0,0,32)}
family={}
for key,r in rd.items():
 f='crystalline_silicate'
 if key in ['limestone','dolostone','chalk','marble']:f='carbonate'
 if key in ['sandstone_quartz','arkose','greywacke','conglomerate','siltstone']:f='cemented_clastic'
 if key in ['mudstone','shale']:f='argillaceous'
 if key in ['rock_salt','gypsum_rock']:f='evaporite'
 if key=='obsidian':f='glass'
 if key in ['scoria','pumice','ignimbrite','tuff_ash']:f='porous_volcanic'
 if r['fabric'] in ['cleaved','foliated']:f='foliated'
 if key=='coal':f='organic'
 if key in ['laterite','gossan']:f='residual'
 if key in ['gravel','sand','silt_clay','till','regolith']:f='unconsolidated'
 family[key]=f
 t,e,k,fr=chains[f]
 if key=='gypsum_rock':e=1500 # 30 GPa for 20 MPa UCS; distinct from salt
 if key=='chalk':e=500 # 4 GPa porous carbonate
 if f!='unconsolidated':
  u=r['mechanical']['ucs_pa']['estimate']['value']
  for name,v,unit in [('tensile_pa',t*u,'Pa'),('young_modulus_pa',e*u,'Pa'),('fracture_toughness_pa_sqrt_m',k*1e6,'Pa sqrt(m)')]:r['mechanical'][name]=dict(status='known_prior',estimate=est(v,unit))
 r['mechanical']['friction_angle_rad']=dict(status='known_prior',estimate=q(math.radians(fr),math.radians(max(5,fr-10)),math.radians(fr+10),'rad','dry family friction prior; saturation, effective stress and joint infill unresolved'))
# F5: resolved zero-aperture fabric; joints are genuine conditional discontinuities only.
fabrics={}; mechanics={}
spacing={'crystalline_silicate':(1.2,8.,.00015),'carbonate':(.6,4.,.0002),'cemented_clastic':(.7,4.,.0002),'argillaceous':(.25,2.,.0001),'evaporite':(2.,6.,.00005),'glass':(.45,2.,.0001),'porous_volcanic':(.35,2.,.0002),'foliated':(.8,5.,.00015),'organic':(.15,1.,.0001),'residual':(.3,1.,.0002)}
for key,r in rd.items():
 f=family[key];ty='none';s=.001;t=k=1.;pen=False
 if r['fabric'] in ['bedded'] or key in ['mudstone','chert']:ty='bedding';s=.12
 if key=='slate':ty='slaty_cleavage';s=.001;t=.25;k=.35;pen=True
 if key in ['phyllite','schist']:ty='schistosity';s=.002 if key=='phyllite' else .005;t=.3;k=.4;pen=True
 if key in ['gneiss','amphibolite','serpentinite','mylonite']:ty='gneissic_banding' if key=='gneiss' else 'foliation';s=.04 if key=='gneiss' else .005;t=.5;k=.6;pen=True
 if key=='obsidian':ty='flow_banding';s=.03
 if key=='basalt_tholeiitic':ty='columnar';s=.6
 if key=='shale':ty='bedding';s=.002;t=.35;k=.45;pen=True
 fabrics[key]=dict(kind=ty,resolved_spacing_m=None if ty=='none' else est(s,'m'),tensile_ratio=q(t,t*.5,min(1,t*1.4)),toughness_ratio=q(k,k*.5,min(1,k*1.4)),aperture_m=0.,penetrative=pen)
 t,e,kic,fr=chains[f]
 if key=='gypsum_rock':e=1500
 if key=='chalk':e=500
 mechanics[key]=dict(family=f,tensile_ucs_ratio=est(t) if t else None,young_ucs_ratio=est(e) if e else None,creep=key in ['rock_salt','gypsum_rock','shale'],creep_conditions='time-dependent response possible; salt pressure/temperature, gypsum stress/moisture, shale clay/water; no creep-rate law supplied' if key in ['rock_salt','gypsum_rock','shale'] else 'creep not represented by this conditional prior; false is not proof of absence')
 r['joint_sets']=[]
 if f=='unconsolidated':continue
 sp,pe,ap=spacing[f]
 if key in ['chert','banded_iron_formation']:sp,pe=.15,1.5
 if key=='basalt_tholeiitic':sp,pe=.6,6.
 desc=[('cross_fabric',90,0),('cross_fabric',90,90)]
 if ty=='bedding':desc+=[('bedding_discontinuity',0,0)]
 if ty=='columnar':desc=[('cooling_column_side',90,a) for a in [0,60,120]]+[('cooling_cross_joint',0,0)]
 for rel,pol,az in desc:
  space=max(s,sp*.4) if rel=='bedding_discontinuity' else sp
  r['joint_sets'].append(dict(relation=rel,spacing_m=q(space,space*.15,space*5,'m'),persistence_m=q(pe,pe*.1,pe*6,'m'),aperture_m=q(ap,0,ap*20,'m'),normal_polar_deg=q(pol,max(0,pol-15),min(180,pol+15),'deg'),normal_azimuth_deg=q(az,max(0,az-15),min(360,az+15),'deg'),reference='RS3',qualification='conditional discontinuity prior by mechanical family; not penetrative fabric and not observed site geometry'))
# Sandstone F2: replace fictitious biotite with matrix/cement/coating inventory.
for key,m in [('sandstone_quartz',dict(quartz=.935,k_feldspar=.015,illite=.025,goethite=.025)),('arkose',dict(quartz=.57,k_feldspar=.30,plagioclase=.05,illite=.055,goethite=.025))]:rd[key]['mode']=mode(m)
def density(r):
 m={x['mineral']:x['mass_fraction']['value'] for x in r['mode']};rho=(1-r['intrinsic_porosity']['value'])/sum(v/md[k]['density_kg_m3']['value'] for k,v in m.items());r['dry_bulk_density_kg_m3']=est(rho,'kg/m3','derived harmonic mass-mode density times solid fraction');r['dry_deposit_density_kg_m3']=est(rho*(1-r['interparticle_porosity']['value']),'kg/m3');r['crystallinity']=q(1-sum(v for k,v in m.items() if k in ['volcanic_glass','felsic_glass_low_absorption','organic_matter']),0,1)
for r in rd.values():density(r)
real=[]
def add(key,name,settings,m=None,params=None,por=None,grain=None,ucs=None,fabric=None):
 r=copy.deepcopy(rd[key]);p={} if params is None else params
 if m is not None:r['mode']=mode(m)
 if por is not None:r['intrinsic_porosity']=q(por,por*.5,min(.95,por*1.5));r['vesicularity']=q(min(r['vesicularity']['value'],por),0,por)
 if grain is not None:
  r['optical_domain_diameter_m']=est(grain,'m')
  for field,v in [('d10_m',grain*.2),('d50_m',grain),('d90_m',grain*4)]:r['grain'][field]=est(v,'m')
 if ucs is not None:
  old_u=r['mechanical']['ucs_pa']['estimate']['value'];ratio=ucs*1e6/old_u
  for n in ['ucs_pa','tensile_pa','young_modulus_pa']:r['mechanical'][n]['estimate']=est(r['mechanical'][n]['estimate']['value']*ratio,'Pa')
 density(r)
 real.append(dict(id=key+'.'+name,lithofacies=key,settings=settings,parameters=p,spec=r,fabric=copy.deepcopy(fabric or fabrics[key]),mechanics=mechanics[key]))
for key in rd:
 p={}
 if key=='shale':p={'organic_mass_fraction':.05}
 if key=='chert':p={'chert_variant':'grey'}
 if key=='ignimbrite':p={'welding_grade':0.}
 if key in ['sandstone_quartz','arkose']:
  p={'sedimentary':dict(maturity=.95 if key=='sandstone_quartz' else .5,cement_mineral='quartz',cement_mass_fraction=.10,grain_coatings=[dict(mineral='goethite',mass_fraction=.025)])}
 add(key,'reference',['declared_reference'],params=p)
for subtype,m in [('i_type',dict(quartz=.28,k_feldspar=.32,plagioclase=.30,biotite=.06,hornblende=.04)),('s_type',dict(quartz=.32,k_feldspar=.32,plagioclase=.24,muscovite=.08,biotite=.04))]:
 for color in ['pink','white']:
  add('granite',subtype+'_'+color,['continental_arc' if subtype=='i_type' else 'continental_collision'],m,dict(granite_type=subtype,k_feldspar_color=color,plagioclase_an_fraction=.25 if subtype=='i_type' else .15))
add('granite','a_type_white',['continental_rift'],dict(quartz=.30,k_feldspar=.48,plagioclase=.15,biotite=.05,magnetite=.02),dict(granite_type='a_type',k_feldspar_color='white',plagioclase_an_fraction=.10))
add('andesite','arc_hornblende',['arc_hydrous'],params=dict(plagioclase_an_fraction=.45))
add('andesite','two_pyroxene',['arc_drier'],dict(plagioclase=.58,augite=.18,orthopyroxene=.14,volcanic_glass=.07,magnetite=.03),dict(plagioclase_an_fraction=.55))
add('dacite','arc_hornblende',['arc_hydrous'],params=dict(plagioclase_an_fraction=.35))
add('dacite','two_pyroxene',['arc_drier'],dict(plagioclase=.50,quartz=.20,augite=.08,orthopyroxene=.07,felsic_glass_low_absorption=.15),dict(plagioclase_an_fraction=.4))
add('basalt_tholeiitic','lip',['large_igneous_province'],params=dict(plagioclase_an_fraction=.6),por=.04)
add('basalt_tholeiitic','island_shield',['ocean_island_shield'],params=dict(plagioclase_an_fraction=.6),por=.05)
add('basalt_alkali','island_cap',['ocean_island_cap'],dict(plagioclase=.35,augite=.28,olivine=.28,magnetite=.04,volcanic_glass=.05),dict(plagioclase_an_fraction=.5),por=.10,grain=.0006)
add('chert','flint',['silica_diagenesis'],dict(quartz=.96,organic_matter=.04),dict(chert_variant='flint'))
add('chert','white',['silica_diagenesis'],dict(quartz=1.),dict(chert_variant='white'))
add('ignimbrite','unwelded',['pyroclastic_flow'],params=dict(welding_grade=0.))
f=copy.deepcopy(fabrics['ignimbrite']);f.update(kind='flow_banding',penetrative=True,resolved_spacing_m=est(.008,'m'),tensile_ratio=q(.6,.3,.9),toughness_ratio=q(.7,.4,1.))
# Same composition: welding densifies and reduces effective optical boundary scattering.
add('ignimbrite','welded',['pyroclastic_flow'],params=dict(welding_grade=1.),por=.08,ucs=90,fabric=f)
real[-1]['spec']['structure_v2']['optical_boundary_efficiency']=q(.025,.005,.1,cond='welding removes ash boundaries; connected glass transport prior, not a fitted tint')
real[-1]['spec']['structure_v2']['kind']='banded';real[-1]['spec']['structure_v2']['wavelength_m']=est(.008,'m')
real[-1]['spec']['joint_sets']=copy.deepcopy(rd['basalt_tholeiitic']['joint_sets'])
add('shale','grey',['marine_clastic'],dict(illite=.47,smectite=.15,quartz=.37,organic_matter=.01),dict(organic_mass_fraction=.01))
add('shale','organic_rich',['anoxic_marine'],dict(illite=.40,smectite=.15,quartz=.30,organic_matter=.15),dict(organic_mass_fraction=.15))
# Inventory group labels require an explicit species split; no silent synonym claim.
ids=['silicic_glass','quartz','k_feldspar','plagioclase','biotite','clinopyroxene','olivine','amphibole','calcite','dolomite','clay','muscovite','halite','gypsum']
groups={'silicic_glass':['felsic_glass_low_absorption','volcanic_glass'],'clinopyroxene':['augite'],'amphibole':['hornblende'],'clay':['illite','smectite','kaolinite']}
mapping=[dict(inventory_id=k,rockspec_ids=groups.get(k,[k]),requires_species_resolution=k in groups,note='group requires body-admitted species split; candidates are supported catalogue species, not exhaustive mineralogy' if k in groups else 'identity mapping; retain source density in conversion receipt') for k in ids]
plumes=[dict(setting=s,weights={k:q(v,0,1,cond='conditional scenario weights, not measured frequency; resolve once per admitted body') for k,v in w.items()}) for s,w in [('large_igneous_province',{'basalt_tholeiitic':.97,'basalt_alkali':.03}),('ocean_island_shield',{'basalt_tholeiitic':.90,'basalt_alkali':.10}),('ocean_island_cap',{'basalt_tholeiitic':.20,'basalt_alkali':.80})]]
alteration=[dict(kind=k,product_mode=mode(m),qualification='conditional replacement assemblage; no automatic reaction, element balance or body admission') for k,m in [
 ('potassic',dict(k_feldspar=.65,biotite=.25,quartz=.10)),
 ('phyllic',dict(muscovite=.55,quartz=.40,pyrite=.05)),
 ('argillic',dict(kaolinite=.60,smectite=.30,illite=.10)),
 ('propylitic',dict(chlorite=.40,epidote=.30,calcite=.30)),
 ('weathering',dict(kaolinite=.50,goethite=.30,quartz=.20))]]
mineralization=[dict(id=k,texture=t,product_mode=mode(m),qualification='illustrative body-admitted fill assemblage only; neither deposit frequency nor recoverable reserve') for k,t,m in [
 ('quartz_sulfide_vein','vein_fill',dict(quartz=.70,pyrite=.10,galena=.10,sphalerite=.10)),
 ('porphyry_sulfides','disseminated',dict(quartz=.90,pyrite=.07,chalcopyrite=.03)),
 ('cassiterite_greisen','vein_fill',dict(quartz=.70,muscovite=.25,cassiterite=.05)),
 ('gold_placer','placer',dict(quartz=.999,native_gold=.001))]]
out=dict(schema_version='pf-rockspec-assemblage/3',core=c,mineral_metadata=meta,realizations=real,inventory_mapping=mapping,plume_priors=plumes,alteration_assemblages=alteration,mineralization_assemblages=mineralization)
(ROOT/'data/catalog_v3.json').write_text(json.dumps(out,indent=2)+'\n')
print(f"v3: {len(rd)} lithofacies, {len(md)} minerals, {len(real)} realizations, {len(mapping)} inventory IDs")
