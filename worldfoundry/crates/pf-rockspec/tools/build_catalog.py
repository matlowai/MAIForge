#!/usr/bin/env python3
"""Rebuild versioned project priors. No downloads and no external packages.
References motivate classifications/relationships, NOT measured endpoints.
"""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
refs = {
 'DHZ': ('Deer, Howie & Zussman, An Introduction to the Rock-Forming Minerals, 3rd ed., 2013', 'Mineral density, hardness, cleavage, optical classes and alteration. Broad project envelopes, not transcribed specimen measurements.'),
 'HOM': ('Anthony et al., Handbook of Mineralogy, Mineral Data Publishing, 1990 onward', 'Mineral species descriptions. Offline bibliographic attribution; numeric endpoints not independently checked in this run.'),
 'LEMAITRE': ('Le Maitre (ed.), Igneous Rocks: A Classification and Glossary of Terms, 2nd ed., 2002', 'Classification context only. Does not prescribe these mass-mode priors.'),
 'BOGGS': ('Boggs, Principles of Sedimentology and Stratigraphy, 5th ed., 2012', 'Sedimentary composition, texture and porosity context.'),
 'YARDLEY': ('Yardley, An Introduction to Metamorphic Petrology, 1989', 'Metamorphic mineral assemblage and fabric context.'),
 'GOODMAN': ('Goodman, Introduction to Rock Mechanics, 2nd ed., 1989', 'Rock mechanical-property scale and variability. Chosen engineering prior endpoints require specimen qualification.'),
 'ISRM': ('ISRM, Suggested Methods for the Quantitative Description of Discontinuities in Rock Masses, 1978', 'Spacing, persistence, aperture and orientation terminology; site joint populations are not implied.'),
 'BIRKELAND': ('Birkeland, Soils and Geomorphology, 3rd ed., 1999', 'Weathering pathways and secondary mineral context; no rates inferred.'),
 'HAPKE': ('Hapke, Theory of Reflectance and Emittance Spectroscopy, 2nd ed., 2012', 'Particulate transport and grain-size dependence; RGB coefficients here are broad project priors, not spectra.'),
 'KM': ('Kubelka & Munk, Ein Beitrag zur Optik der Farbanstriche, 1931', 'Two-flux diffuse half-space R=(1+K/S)-sqrt((1+K/S)^2-1). Intimate mixture used as a reduced model, not a crystalline BRDF solution.'),
 'LD': ('Lekner & Dorf, Why some things are darker when wet, Applied Optics 27, 1278-1280, 1988', 'Escape-cone internal reflection mechanism. This implementation is a declared reduced multiple-return model, not a reproduction of all paper terms.'),
 'FRESNEL': ('Born & Wolf, Principles of Optics, 7th ed., 1999', 'Normal-incidence dielectric/complex-index Fresnel equations.'),
 'PROJECT': ('CX-MAT-ROCKSPEC v1, operator brief 2026-10-25 and README.md', 'Explicit project prior or model closure. No measured confidence interval. Covariance groups denote constraints, not a fitted covariance matrix.'),
 'BRIEF': ('design CX-MAT-ROCKSPEC brief, 2026-10-25', 'Supplied sanity examples basalt 0.05-0.10, granite 0.25-0.35, white marble >=0.6. Not independently verified literature measurements.')
}
def q(v, lo, hi, unit, ref='PROJECT', cond='ambient dry representative material; prior not measurement', group='material'):
 return dict(value=v, interval=[lo,hi], unit=unit, interval_kind='plausible_prior_support', distribution='unspecified_joint_prior', reference=ref, evidence_kind='project_prior', support_scale='specimen; not a site estimate', conditions=cond, covariance_group=group)
def frac(v, ref='PROJECT', group='mode'):
 return q(v,max(0,v*.55),min(1,v*1.45), '1',ref, 'dry solid mass fraction; jointly constrained to unit sum',group)
def interval(v, unit, ref, low=.6, high=1.5, **kw): return q(v,v*low,v*high,unit,ref,**kw)
def known(e): return {'status':'known_prior','estimate':e}
def unknown(s): return {'status':'unknown','reason':s}
# id, rho kg/m3, hardness, n, cleavage directions, angles deg, quality,
# fracture habit, lustre, diffuse linear RGB at 1 mm particulate support, products.
rows = [
('quartz',2650,7,1.55,0,[],'none','conchoidal','vitreous',[.78,.78,.77],[]),
('plagioclase',2700,6,1.56,2,[86,94],'good','uneven','vitreous',[.58,.59,.60],['kaolinite','smectite']),
('k_feldspar',2570,6,1.52,2,[90],'good','uneven','vitreous',[.64,.47,.39],['kaolinite','illite']),
('biotite',3100,2.75,1.61,1,[],'perfect','uneven','pearly',[.035,.026,.019],['chlorite','vermiculite','goethite']),
('muscovite',2830,2.5,1.59,1,[],'perfect','uneven','pearly',[.68,.65,.57],['illite','kaolinite']),
('hornblende',3250,5.5,1.68,2,[56,124],'good','uneven','vitreous',[.025,.035,.029],['chlorite','goethite']),
('augite',3400,5.5,1.71,2,[87,93],'good','uneven','vitreous',[.025,.032,.024],['chlorite','smectite','goethite']),
('olivine',3350,6.75,1.70,2,[90],'poor','conchoidal','vitreous',[.12,.16,.065],['serpentine','iddingsite']),
('orthopyroxene',3350,5.5,1.70,2,[88,92],'good','uneven','vitreous',[.065,.050,.035],['serpentine','goethite']),
('calcite',2710,3,1.57,3,[75,105],'perfect','uneven','vitreous',[.78,.77,.73],['dissolved_carbonate']),
('dolomite',2860,3.75,1.59,3,[74,106],'perfect','uneven','vitreous',[.67,.64,.56],['calcite','dissolved_carbonate']),
('halite',2160,2.5,1.54,3,[90],'perfect','conchoidal','vitreous',[.77,.78,.77],['dissolved_salt']),
('gypsum',2320,2,1.52,3,[66,114],'one_perfect_two_weaker','uneven','pearly',[.76,.74,.68],['dissolved_sulfate']),
('kaolinite',2600,2.25,1.56,1,[],'perfect_microcrystal','earthy','earthy',[.71,.67,.58],[]),
('smectite',2350,1.5,1.52,1,[],'perfect_microcrystal','earthy','earthy',[.36,.32,.23],[]),
('illite',2750,1.5,1.56,1,[],'perfect_microcrystal','earthy','earthy',[.40,.36,.28],['kaolinite']),
('chlorite',2850,2.5,1.60,1,[],'perfect','uneven','pearly',[.12,.19,.11],['smectite','goethite']),
('serpentine',2550,3.5,1.55,1,[],'variable_polymorph','splintery','waxy',[.15,.24,.12],['smectite','goethite']),
('garnet',3900,7,1.80,0,[],'none','conchoidal','vitreous',[.19,.045,.035],['goethite','kaolinite']),
('magnetite',5180,5.75,2.4,0,[],'none','uneven','sub_metallic',[.012,.012,.013],['hematite','goethite']),
('ilmenite',4720,5.5,2.5,0,[],'none','conchoidal','sub_metallic',[.018,.017,.018],['leucoxene','hematite']),
('hematite',5260,6,2.9,0,[],'none','uneven','sub_metallic',[.22,.042,.025],['goethite']),
('goethite',4280,5.25,2.3,1,[],'perfect_microcrystal','uneven','earthy',[.32,.14,.038],['hematite']),
('pyrite',5010,6.25,2.9,3,[90],'indistinct','conchoidal','metallic',[.025,.020,.009],['goethite','sulfate']),
('chalcopyrite',4200,3.75,2.7,2,[90],'poor','uneven','metallic',[.025,.016,.005],['malachite','goethite']),
('volcanic_glass',2400,5.5,1.50,0,[],'none','conchoidal','glassy',[.014,.012,.010],['hydrated_glass','smectite']),
('felsic_glass_low_absorption',2350,5.5,1.49,0,[],'none','conchoidal','glassy',[.69,.71,.72],['hydrated_glass','smectite']),
('organic_matter',1350,2.5,1.70,0,[],'none','conchoidal','dull_to_vitreous',[.020,.018,.016],['oxidized_organic_matter']),
('epidote',3450,6.5,1.75,1,[],'perfect','uneven','vitreous',[.15,.21,.06],['clay_residue','goethite']),
('talc',2780,1,1.56,1,[],'perfect','uneven','pearly',[.61,.66,.56],[]),
('apatite',3200,5,1.64,1,[],'poor','conchoidal','vitreous',[.42,.49,.34],['dissolved_phosphate'])]
minerals=[]
for key,rho,h,n,d,a,c,f,l,rgb,products in rows:
 metals={'magnetite':([2.4]*3,[2.0]*3), 'ilmenite':([2.5]*3,[1.7]*3), 'hematite':([2.9]*3,[.65,.8,1.0]), 'pyrite':([2.9,3.0,3.2],[3.4,2.8,1.9]), 'chalcopyrite':([2.7,2.8,3.0],[3.6,2.6,1.5])}
 minerals.append(dict(id=key,density_kg_m3=interval(rho,'kg/m3','DHZ',.94,1.08), mohs=q(h,max(1,h-.5),min(10,h+.5),'Mohs','HOM'),
  cleavage=dict(directions=q(d,d,d,'count','DHZ',cond='representative species; group/polymorph variation described by quality'),angles_deg=[q(x,x-3,x+3,'deg','DHZ') for x in a],quality=c,reference='DHZ'),
  fracture_habit=f,lustre=l,refractive_index=q(n,max(1,n-.08),n+.12,'1','DHZ',cond='effective visible n; anisotropy/composition unresolved; opaque species require k'),
  extinction_k=None if key not in metals else [q(x,max(0,x*.5),x*1.5,'1','PROJECT',cond='opaque optical surrogate, not measured n/k spectrum') for x in metals[key][1]],
  diffuse_reflectance=[q(x,max(0.001,x*.45),min(.95,x*1.6),'linear_reflectance','HAPKE',cond='clean optically thick granular aggregate, reference diameter 1 mm; not polished crystal face reflectance') for x in rgb],
  optical_reference_diameter_m=q(.001,.001,.001,'m','PROJECT'), internal_scattering_length_m=q(.00005,.000005,.005,'m','PROJECT',cond='effective microcrack/internal scattering length; no calibrated spectrum'),
  alteration_products=products, plagioclase_an_fraction=q(.5,0,1,'1','DHZ',cond='An molar fraction 0 albite to 1 anorthite; resolve n/rho jointly for specific composition') if key=='plagioclase' else None))
# key | modal percentages (solid mass, NOT QAPF volume) | median grain mm |
# porosity | UCS MPa | fabric | weathering family | products
rockrows='''basalt_tholeiitic|plagioclase:48,augite:35,olivine:5,magnetite:5,volcanic_glass:7|0.2|0.05|180|columnar|mafic_oxidation|goethite,smectite
basalt_alkali|plagioclase:40,augite:32,olivine:18,magnetite:5,volcanic_glass:5|0.3|0.06|160|massive|mafic_oxidation|goethite,smectite
andesite|plagioclase:58,hornblende:20,augite:12,volcanic_glass:7,magnetite:3|0.4|0.08|130|massive|silicate_hydrolysis|smectite,goethite
dacite|plagioclase:50,quartz:20,hornblende:10,biotite:5,felsic_glass_low_absorption:15|0.3|0.08|100|massive|silicate_hydrolysis|kaolinite,goethite
rhyolite|quartz:28,k_feldspar:40,plagioclase:15,biotite:2,felsic_glass_low_absorption:15|0.1|0.10|100|massive|silicate_hydrolysis|kaolinite
obsidian|volcanic_glass:98,magnetite:2|0.5|0.01|120|massive|glass_hydration|hydrated_glass,perlite
scoria|plagioclase:35,augite:25,olivine:10,volcanic_glass:25,magnetite:5|0.3|0.55|15|vesicular|mafic_oxidation|goethite,smectite
pumice|felsic_glass_low_absorption:85,quartz:5,k_feldspar:10|0.03|0.75|3|vesicular|glass_hydration|hydrated_glass,smectite
ignimbrite|felsic_glass_low_absorption:55,quartz:15,k_feldspar:20,plagioclase:10|0.5|0.25|35|bedded|glass_hydration|smectite,kaolinite
tuff_ash|felsic_glass_low_absorption:65,plagioclase:20,quartz:10,biotite:5|0.06|0.35|15|clastic|glass_hydration|smectite,kaolinite
granite|quartz:30,k_feldspar:35,plagioclase:30,biotite:5|2|0.01|180|massive|grus_disaggregation|kaolinite,goethite
granodiorite|quartz:25,k_feldspar:15,plagioclase:48,biotite:7,hornblende:5|2|0.015|160|massive|grus_disaggregation|kaolinite,goethite
diorite|plagioclase:65,hornblende:25,biotite:5,quartz:5|2|0.015|170|massive|silicate_hydrolysis|kaolinite,goethite
gabbro|plagioclase:50,augite:40,olivine:7,magnetite:3|3|0.01|200|massive|mafic_oxidation|goethite,smectite
peridotite|olivine:75,orthopyroxene:20,magnetite:5|3|0.01|200|massive|serpentinization|serpentine,goethite
sandstone_quartz|quartz:94,k_feldspar:3,illite:3|0.3|0.18|60|bedded|granular_disaggregation|quartz,goethite
arkose|quartz:60,k_feldspar:30,plagioclase:5,illite:5|0.6|0.18|50|bedded|silicate_hydrolysis|kaolinite,goethite
greywacke|quartz:40,plagioclase:25,chlorite:15,illite:15,biotite:5|0.3|0.08|100|bedded|silicate_hydrolysis|illite,goethite
conglomerate|quartz:55,plagioclase:20,k_feldspar:15,illite:10|20|0.12|70|clastic|granular_disaggregation|quartz,kaolinite
siltstone|quartz:55,illite:30,plagioclase:15|0.025|0.12|60|bedded|slaking|illite,goethite
mudstone|illite:45,smectite:20,quartz:30,organic_matter:5|0.005|0.15|30|massive|slaking|illite,smectite
shale|illite:45,smectite:15,quartz:35,organic_matter:5|0.004|0.12|40|bedded|slaking|illite,goethite
limestone|calcite:95,illite:3,quartz:2|0.2|0.08|80|bedded|carbonate_solution|calcite,clay_residue
dolostone|dolomite:95,calcite:3,quartz:2|0.2|0.08|110|bedded|carbonate_solution|dolomite,clay_residue
chalk|calcite:98,illite:2|0.003|0.40|8|bedded|carbonate_solution|calcite,clay_residue
chert|quartz:98,hematite:2|0.002|0.01|200|massive|silica_residual|quartz
rock_salt|halite:98,gypsum:2|3|0.01|25|massive|evaporite_solution|halite
gypsum_rock|gypsum:95,calcite:5|0.5|0.05|20|bedded|evaporite_solution|gypsum
coal|organic_matter:90,illite:8,pyrite:2|0.1|0.08|20|bedded|organic_oxidation|oxidized_organic_matter,goethite
slate|illite:45,quartz:30,chlorite:20,organic_matter:5|0.02|0.015|120|cleaved|silicate_hydrolysis|illite,goethite
phyllite|muscovite:35,quartz:35,chlorite:25,biotite:5|0.08|0.02|100|foliated|silicate_hydrolysis|illite,goethite
schist|muscovite:35,quartz:30,biotite:20,garnet:10,plagioclase:5|1.5|0.02|100|foliated|silicate_hydrolysis|illite,goethite
gneiss|quartz:30,k_feldspar:25,plagioclase:25,biotite:15,hornblende:5|2|0.01|170|foliated|grus_disaggregation|kaolinite,goethite
marble|calcite:99,dolomite:1|1|0.005|90|massive|carbonate_solution|calcite
quartzite|quartz:98,muscovite:2|0.2|0.005|220|massive|silica_residual|quartz
serpentinite|serpentine:85,magnetite:8,talc:5,olivine:2|0.1|0.03|70|foliated|mafic_oxidation|goethite,smectite
amphibolite|hornblende:60,plagioclase:35,quartz:5|1|0.01|190|foliated|silicate_hydrolysis|chlorite,goethite
hornfels|quartz:40,plagioclase:25,biotite:20,cordierite:15|0.1|0.01|200|massive|silicate_hydrolysis|kaolinite,goethite
gravel|quartz:55,plagioclase:20,k_feldspar:15,illite:10|20|0.30|0|clastic|granular_disaggregation|quartz,kaolinite
sand|quartz:85,k_feldspar:10,illite:5|0.4|0.36|0|clastic|granular_disaggregation|quartz,kaolinite
silt_clay|illite:40,smectite:25,quartz:30,kaolinite:5|0.006|0.45|0|clastic|slaking|illite,smectite
till|quartz:40,plagioclase:20,k_feldspar:15,illite:15,chlorite:10|3|0.25|0|clastic|silicate_hydrolysis|kaolinite,goethite
regolith|quartz:40,kaolinite:30,illite:15,goethite:10,plagioclase:5|0.5|0.40|0|clastic|residual_weathering|kaolinite,goethite'''
# Hornfels needs a real accessory instead of an invented feldspar substitution.
rows_extra=('cordierite',2600,7.25,1.55,1,[],'poor','conchoidal','vitreous',[.32,.36,.40],['illite','chlorite'])
key,rho,h,n,d,a,c,f,l,rgb,products=rows_extra
m=json.loads(json.dumps(minerals[0]));m.update(id=key,density_kg_m3=interval(rho,'kg/m3','DHZ',.94,1.08),mohs=q(h,7,7.5,'Mohs','HOM'),refractive_index=q(n,1.52,1.58,'1','DHZ'),cleavage=dict(directions=q(d,d,d,'count','DHZ'),angles_deg=[],quality=c,reference='DHZ'),fracture_habit=f,lustre=l,alteration_products=products,diffuse_reflectance=[q(x,x*.45,min(.95,x*1.6),'linear_reflectance','HAPKE',cond='clean optically thick granular aggregate at 1 mm; optical prior') for x in rgb]);minerals.append(m)
md={m['id']:m for m in minerals}
rocks=[]
for i,line in enumerate(rockrows.splitlines()):
 key,mode,grain,por,ucs,fabric,family,prod=line.split('|'); grain=float(grain)/1000;por=float(por);ucs=float(ucs)*1e6
 mode={k:float(v)/100 for k,v in (x.split(':') for x in mode.split(','))}
 ref='LEMAITRE' if i<15 else 'BOGGS' if i<29 or i>=38 else 'YARDLEY'
 loose=i>=38; ves=por*.95 if key in ('scoria','pumice') else por*.2 if i<10 else 0
 intrinsic=0 if loose else por; deposit=por if loose else 0
 density=(1-intrinsic)/sum(v/md[k]['density_kg_m3']['value'] for k,v in mode.items())
 mech={k:known(interval(v,u,'GOODMAN',.25,2.5,cond='intact dry specimen engineering prior; anisotropy/confinement not resolved')) for k,v,u in [('ucs_pa',ucs,'Pa'),('tensile_pa',ucs*.065,'Pa'),('young_modulus_pa',ucs*300,'Pa'),('fracture_toughness_pa_sqrt_m',max(.1e6,ucs*.008),'Pa sqrt(m)')]}
 if loose:
  mech={k:unknown('unconsolidated aggregate; intact-rock parameter not applicable; constitutive effective-stress model required') for k in mech}
 mech['poisson_ratio']=known(q(.25,.15,.4,'1','GOODMAN')) if not loose else unknown('effective stress and packing dependent; not an intact elastic solid')
 mech['friction_angle_rad']=known(q(.60,.35,.85,'rad','GOODMAN',cond='dry friction prior; no cohesion inferred; effective stress dependent'))
 joints=[]
 if not loose:
  desc=[('fabric_parallel',0,0)] if fabric in ('bedded','foliated','cleaved') else [('cross_fabric',90,0)]
  desc += [('cross_fabric',90,90)]
  if fabric=='columnar': desc=[('cooling_column_side',90,0),('cooling_column_side',90,60),('cooling_column_side',90,120),('cooling_cross_joint',0,0)]
  for relation,dip,az in desc:
   space=.03 if fabric=='cleaved' else .15 if fabric=='foliated' else .5
   joints.append(dict(relation=relation,spacing_m=q(space,space*.1,space*10,'m','ISRM'),persistence_m=q(5,.2,50,'m','ISRM'),aperture_m=q(.0003,0,.01,'m','ISRM'),normal_polar_deg=q(dip,max(0,dip-20),min(180,dip+20),'deg','PROJECT'),normal_azimuth_deg=q(az,max(0,az-20),min(360,az+20),'deg','PROJECT'),reference='ISRM',qualification='conditional set prior, not a generated fracture or observed site; normal in fabric coordinates'))
 rocks.append(dict(id=key,mode=[dict(mineral=k,mass_fraction=frac(v,ref,key+'_closed_mode')) for k,v in mode.items()],
 grain=dict(d10_m=q(grain*.2,grain*.02,grain*.5,'m',ref),d50_m=q(grain,grain*.5,grain*2,'m',ref),d90_m=q(grain*4,grain*2,grain*20,'m',ref),distribution='lognormal_descriptor_prior; percentiles not independently sampled',aspect_ratio=q(3 if fabric in ('foliated','cleaved') else 1.5,1,10,'1',ref),basis='clast diameter; internal crystal size unresolved' if key in ('conglomerate','gravel','till') else 'glass defect/inclusion spacing, not crystal diameter' if key=='obsidian' else 'mineral/crystal/particle diameter'),
 optical_domain_diameter_m=q(.0005 if key in ('conglomerate','gravel','till') else grain,.00001 if key in ('conglomerate','gravel','till') else grain*.2,.01 if key in ('conglomerate','gravel','till') else grain*4,'m','PROJECT',cond='mineral/internal optical scattering domain; separate from clast size in gravel, till, conglomerate'),
 crystallinity=q(max(0,1-mode.get('volcanic_glass',0)-mode.get('felsic_glass_low_absorption',0)-mode.get('organic_matter',0)),0,1,'1',ref,cond='crystalline fraction of solids; glass and organic matter treated disordered'),vesicularity=q(ves,0,min(.95,ves*1.3),'1',ref,cond='subset of intrinsic porosity; unresolved voids only'),intrinsic_porosity=q(intrinsic, max(0,intrinsic*.5),min(.95,max(intrinsic*1.4,.02 if not loose else 0)),'1',ref),interparticle_porosity=q(deposit,deposit*.7,min(.75,deposit*1.3),'1',ref),
 dry_bulk_density_kg_m3=q(density,density*.8,density*1.15,'kg/m3','PROJECT',cond='derived nominal (1-intrinsic_porosity)/sum(w/rho); excludes deposit interparticle pores and resolved cavities',group=key+'_mode_porosity_density'),
 dry_deposit_density_kg_m3=q(density*(1-deposit),density*(1-deposit)*.7,density*(1-deposit)*1.25,'kg/m3','PROJECT',cond='intrinsic dry bulk density times (1-interparticle_porosity); no resolved cavities',group=key+'_mode_porosity_density'),
 mechanical=mech,fabric=fabric,joint_sets=joints,weathering_family=family,weathering_products=prod.split(','),reference=ref,
 thermal_response=unknown('thermal law/calibration outside this lane; do not substitute generic conductivity'),weathering_rate=unknown('pathways only; time-dependent rates outside scope')))
model={
 'water_ior':q(1.333,1.33,1.34,'1','FRESNEL',cond='visible liquid water near ambient conditions'),
 'roughness_by_lustre': {k:q(v,max(.02,v-.2),min(1,v+.2),'GGX_perceptual_roughness','PROJECT',cond='uncalibrated fracture-surface microfacet prior') for k,v in dict(vitreous=.48,pearly=.35,sub_metallic=.42,metallic=.3,glassy=.14,earthy=.9,waxy=.5,dull_to_vitreous=.5).items()},
 'cleavage_slope':q(.18,.03,.4,'radian_rms','PROJECT'), 'fracture_slope':q(.38,.1,.8,'radian_rms','PROJECT'),
 'glass_slope':q(.04,.005,.2,'radian_rms','PROJECT',cond='smooth conchoidal glass faces; excludes resolved curved fracture geometry'),
 'fabric_slope_reduction':q(.5,0,.8,'1','PROJECT'),
 'sparkle_fraction':q(.12,0,.3,'area_fraction','PROJECT',cond='upper cleavage-face coverage prior times eligible mineral volume fraction; redistribution of specular energy'),
 'film_smoothing_length_m':q(.00005,.000005,.0005,'m','PROJECT'),
 'rind_extinction_m_inv':q(10000,1000,100000,'1/m','PROJECT',cond='effective two-pass attenuation through optically rough coating; uncalibrated'),
}
catalog=dict(schema_version='pf-rockspec/1',numerical_model='km-internal-scattering+escape-cone/1', references={k:dict(citation=v[0],qualification=v[1],verification='offline_bibliographic_reference_not_source_extraction') for k,v in refs.items()},minerals=minerals,rocks=rocks,model=model,
 sanity_targets=[dict(rock=k,lower=lo,upper=hi,reference='BRIEF',qualification='operator-supplied approximate sanity range; not independently verified literature') for k,lo,hi in [('basalt_tholeiitic',.05,.10),('granite',.25,.35),('marble',.6,.9)]])
(ROOT/'data/catalog_v1.json').write_text(json.dumps(catalog,indent=2)+'\n')
(ROOT/'data/shared_keys_v1.json').write_text(json.dumps([r['id'] for r in rocks],indent=2)+'\n')
print(f'{len(minerals)} minerals, {len(rocks)} rocks')
