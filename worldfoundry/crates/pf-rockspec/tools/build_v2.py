#!/usr/bin/env python3
"""Generate v2 conditional priors from immutable v1. Standard library, offline.
No RGB tint fitting. Changed modes and connected-boundary efficiencies are
explicit new hypotheses, not measured specimen composition.
"""
import copy,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'data/catalog_v1.json').read_text())
def q(v,lo,hi,unit,cond):
 return dict(value=v,interval=[lo,hi],unit=unit,interval_kind='plausible_prior_support',distribution='unspecified_joint_prior',reference='RS2',evidence_kind='project_prior',support_scale='conditional specimen, not site measurement',conditions=cond,covariance_group='v2_structure_transport')
c['schema_version']='pf-rockspec/2';c['numerical_model']='connected-km+cavity-escape+material-structure/2'
c['references']['RS2']=dict(citation='CX-MAT-ROCKSPEC-2 model and conditional priors, 2026-10-25; MODEL_V2.md',qualification='Explicit uncalibrated physical closures and alternative assemblages; not source-extracted endpoints',verification='project_model')
# Mn oxide is an optical constituent in explicitly supplied arid varnish only.
mn=copy.deepcopy(next(m for m in c['minerals'] if m['id']=='hematite'));mn['id']='pyrolusite'
mn['density_kg_m3']=q(5000,4700,5200,'kg/m3','MnO2 conditional mineral prior')
mn['refractive_index']=q(2.4,2.,3.,'1','effective visible n; not measured dispersion')
mn['extinction_k']=[q(1.1,.5,2.,'1','broad effective extinction prior') for _ in range(3)]
mn['diffuse_reflectance']=[q(v,.005,.09,'linear_reflectance','dark Mn oxide aggregate prior') for v in [.026,.023,.020]]
mn['mohs']=q(6,2,6.5,'Mohs','crystalline versus earthy habit variability')
mn['cleavage']=dict(directions=q(1,1,1,'count','representative pyrolusite cleavage'),angles_deg=[],quality='perfect',reference='HOM')
mn['alteration_products']=[];mn['lustre']='sub_metallic';c['minerals'].append(mn)
md={m['id']:m for m in c['minerals']}
# Replace part of host mass with declared minor phases. No stock is created by evaluation.
additions={
 'rhyolite':{'magnetite':.01},
 'sandstone_quartz':{'goethite':.025,'biotite':.025,'organic_matter':.004},
 'arkose':{'goethite':.025,'biotite':.04,'organic_matter':.004},
 'greywacke':{'organic_matter':.025,'goethite':.02},
 'conglomerate':{'biotite':.055,'goethite':.025,'organic_matter':.006},
 'siltstone':{'organic_matter':.035,'goethite':.03},
 'limestone':{'organic_matter':.006,'goethite':.009},
 'dolostone':{'organic_matter':.008,'goethite':.01},
 'chert':{'organic_matter':.035,'goethite':.015},
 'phyllite':{'organic_matter':.02},
 'quartzite':{'organic_matter':.006,'goethite':.02},
 'gravel':{'biotite':.055,'goethite':.025,'organic_matter':.006},
 'sand':{'biotite':.025,'goethite':.025,'organic_matter':.004},
 'silt_clay':{'organic_matter':.035,'goethite':.03},
 'till':{'organic_matter':.02,'biotite':.025,'goethite':.025},
 'regolith':{'organic_matter':.02},
}
efficiencies={'shale':.015,'mudstone':.02,'siltstone':.08,'silt_clay':.06,'chert':.008,'slate':.08,'phyllite':.3,'pumice':.15,'chalk':.12,'quartzite':.35}
for r in c['rocks']:
 key=r['id'];old={m['mineral']:m['mass_fraction']['value'] for m in r['mode']}
 extra=additions.get(key,{})
 mode={k:v*(1-sum(extra.values())) for k,v in old.items()}
 for k,v in extra.items():mode[k]=mode.get(k,0)+v
 if key=='gneiss': mode=dict(quartz=.34,k_feldspar=.29,plagioclase=.25,biotite=.08,hornblende=.04)
 r['mode']=[dict(mineral=k,mass_fraction=q(v,max(0,v*.4),min(1,v*1.6),'1','conditional v2 dry mass mode; closed jointly; not field-admitted')) for k,v in mode.items()]
 rho=(1-r['intrinsic_porosity']['value'])/sum(v/md[k]['density_kg_m3']['value'] for k,v in mode.items())
 r['dry_bulk_density_kg_m3']=q(rho,rho*.8,rho*1.2,'kg/m3','derived from v2 modes and intrinsic porosity')
 dep=rho*(1-r['interparticle_porosity']['value']);r['dry_deposit_density_kg_m3']=q(dep,dep*.8,dep*1.2,'kg/m3','derived with separate packing porosity')
 crystal=1-sum(v for k,v in mode.items() if k in ['organic_matter','volcanic_glass','felsic_glass_low_absorption'])
 r['crystallinity']=q(crystal,0,1,'1','crystalline solid fraction from declared mode')
 kind='crystalline';period=r['grain']['d50_m']['value']
 if key in ('gneiss',):kind='banded';period=.06
 elif key in ('shale','siltstone','slate','phyllite'):kind='laminated';period=.003 if key in ('shale','slate') else .006
 elif key in ('conglomerate','gravel','till'):kind='clastic';period=r['grain']['d50_m']['value']
 elif key in ('andesite','dacite'):kind='porphyritic';period=.006;r['grain']['d90_m']=q(.006,.002,.015,'m','conditional phenocryst upper grain population')
 elif r['fabric']=='vesicular':kind='vesicular';period=.004 if key=='pumice' else .008
 elif key=='obsidian':kind='glass'
 elif key=='schist':kind='micaceous'
 if key=='till':r['grain']['d90_m']=q(.08,.02,.5,'m','conditional poorly sorted coarse clast population');period=.02
 eff=efficiencies.get(key,1.)
 r['structure_v2']=dict(kind=kind,optical_boundary_efficiency=q(eff,max(.001,eff*.2),1.,'1','fraction of grain-boundary scattering retained by connected/intergrown or aggregated domains; d_transport=d_optical/efficiency; uncalibrated'),wavelength_m=q(period,period*.3,period*3,'m','band/lamina/vesicle spacing or coarse population support; actual physical texture scale'),contrast=q(.8,0,1,'1','bounded mean-preserving mineral heterogeneity closure, not independent colours'))
c['sanity_targets'].append(dict(rock='obsidian',lower=.008,upper=.03,reference='RS2',qualification='arc-1 dark obsidian preservation envelope, not literature'))
(ROOT/'data/catalog_v2.json').write_text(json.dumps(c,indent=2)+'\n')
print('v2: 43 rocks, 33 constituents; v1 untouched')
