"""USGS-only optical calibration. No RELAB reader or dependency.
Usage: python3 -I build_v5.py COVERAGE_DIR USGS_CVASD_DIR LANE_DIR OUTPUT_JSON
Inputs are read-only text. All fitted observations retain file SHA256 and record.
"""
import sys, json, math, re, hashlib, statistics
from pathlib import Path

def g(x,mu,a,b): return math.exp(-.5*((x-mu)/(a if x<mu else b))**2)
def cmf(x):
    return (1.056*g(x,599.8,37.9,31)-.065*g(x,501.1,20.4,26.2)+.362*g(x,442,16,26.7),
            .821*g(x,568.8,46.9,40.5)+.286*g(x,530.9,16.3,31.1),
            1.217*g(x,437,11.8,36)+.681*g(x,459,26,13.8))
def mat(m,v):return [sum(a*b for a,b in zip(row,v)) for row in m]
B=[[.8951,.2664,-.1614],[-.7502,1.7135,.0367],[.0389,-.0685,1.0296]]
BI=[[.9869929,-.1470543,.1599627],[.4323053,.5183603,.0492912],[-.0085287,.0400428,.9684867]]
RGB=[[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]]
CMF=[cmf(x) for x in range(380,781)]
WHITE=[sum(c[j] for c in CMF) for j in range(3)]
SRC=mat(B,[v/WHITE[1] for v in WHITE]);DST=mat(B,[.95047,1,1.08883])
WEIGHTS=[]
for c in CMF:
    lms=mat(B,[v/WHITE[1] for v in c]);WEIGHTS.append(mat(RGB,mat(BI,[lms[j]*DST[j]/SRC[j] for j in range(3)])))
METHOD='USGS splib07b cvASD; 380-780 nm 1 nm rectangular quadrature; illuminant E; Wyman/Sloan/Shirley 2013 analytic CIE1931 2deg CMF; Bradford E->D65; linear sRGB. RGB reduced transport approximation, not monochromatic coefficients or calibrated BRDF.'

def stats(rows):
    if not rows:return None
    xs=[r['rgb'] for r in rows];n=len(xs);mean=[sum(x[j] for x in xs)/n for j in range(3)]
    covariance=[[sum((x[i]-mean[i])*(x[j]-mean[j]) for x in xs)/n for j in range(3)] for i in range(3)]
    return dict(n_spectra=n,mean=mean,interval=[[min(x[j] for x in xs),max(x[j] for x in xs)] for j in range(3)],covariance=covariance,interval_kind='observed_min_max_not_confidence',spread_kind='population_sd_of_specimen_fraction_spread',sd=[math.sqrt(covariance[j][j]) for j in range(3)],hue_policy='joint RGB observation draws retain cross-channel covariance; no independent hue noise')

def grain(title):
    t=title.lower();m=re.search(r'(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*u(?:m)?\b',t)
    if m:
        a,b=map(float,m.groups());return [a*1e-6,b*1e-6],(a+b)*.5e-6,'title_range_midpoint_assumption'
    m=re.search(r'(?:<|lt)\s*(\d+(?:\.\d+)?)\s*u(?:m)?\b',t)
    if m:return [0,float(m[1])*1e-6],float(m[1])*.5e-6,'title_upper_bound_half_assumption'
    m=re.search(r'\b(\d+(?:\.\d+)?)\s*u(?:m)?\b',t)
    if m:return [float(m[1])*1e-6]*2,float(m[1])*1e-6,'title_nominal'
    return None,None,'unknown_keep_catalogue_reference_diameter_prior'

def main():
    covdir,root,lane,out=map(Path,sys.argv[1:5]);lane.mkdir(parents=True,exist_ok=True)
    cov=json.load(open(covdir/'coverage.json'));idx=json.load(open(covdir/'splib07_index.json'))
    wlfile=next(root.glob('*Wavelengths*'));wl=[float(s) for s in wlfile.read_text().splitlines()[1:] if s.strip()]
    vis=[min(range(len(wl)),key=lambda i:abs(wl[i]*1000-nm)) for nm in range(380,781)]
    assert max(abs(wl[i]*1000-nm) for nm,i in zip(range(380,781),vis))<.01
    rejects=[];cache={}
    def read(r):
        if r['record'] in cache:return cache[r['record']]
        why=None
        if r['spectrometer'].startswith('NIC'):why='no_visible_support'
        elif r['measurement_type']!='AREF':why='relative_reflectance'
        elif r['purity_code'] not in ['a','b','c']:why='purity_unknown_or_significant_contamination'
        p=root/r['file_cvASD'] if r.get('file_cvASD') else None
        if not why and (p is None or not p.is_file()):why='no_cvASD_file'
        if not why:
            vals=[float(x) for x in p.read_text(encoding='latin-1').splitlines()[1:] if x.strip()]
            assert len(vals)==len(wl)
            v=[vals[i] for i in vis]
            if any(not math.isfinite(x) or x<0 or x>1 for x in v):why='invalid_missing_or_nonpassive_visible_value_no_gap_fill'
        if not why:
            rgb=[sum(x*w[j] for x,w in zip(v,WEIGHTS)) for j in range(3)]
            if any(x<=0 or x>=1 for x in rgb):why='out_of_gamut_linear_RGB_no_clipping_for_fit'
        if why:
            rejects.append(dict(record=r['record'],title=r['title'],reason=why));cache[r['record']]=None;return None
        interval,d,assumption=grain(r['title'])
        result=dict(record=r['record'],title=r['title'],sample_id=r['sample_id'],file=r['file_cvASD'],sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rgb=rgb,Y=sum(x*c[1] for x,c in zip(v,CMF))/WHITE[1],grain_m=d,grain_interval_m=interval,grain_basis=assumption,purity=r['purity_code'],form='particulate' if d else 'unspecified')
        if r['record']==11848:result.update(form='solid_fresh_surface',condition='fresh',metadata='https://pubs.usgs.gov/of/2003/ofr-03-395/DESCRIPT/S/basalt_fresh_br93-46b.html')
        if r['record']==12770:result.update(form='solid_weathered_surface',condition='weathered_or_coated',metadata='https://pubs.usgs.gov/of/2003/ofr-03-395/DESCRIPT/S/limestone_cu02-11a.html')
        cache[r['record']]=result;return result
    minerals={};missing=[]
    for m in cov['minerals']:
        groups={}
        for r in sorted(m['usgs_single'],key=lambda r:(0 if r['spectrometer']=='BECK' else 1,r['record'])):
            if r['chapter']!='ChapterM_Minerals' or '+' in r['title'] or '/' in r['title'] or 'Thin_Film' in r['title']:continue
            a=read(r)
            if a:groups.setdefault(r['sample_id'],a)
        rows=list(groups.values())
        if not rows:missing.append(m['id'])
        minerals[m['id']]=dict(observations=rows,statistics=stats(rows),qualification='empirical_reflectance_factor_including_a_b_c_purity_spread' if rows else 'no_usable_USGS_keep_resolved_v3_optical_prior',grain_qualification='Unknown-size observations are not size-calibration points; no HS fraction codes decoded.')
    rocks={}
    for r in cov['rocks']:
        rows=[]
        for e in r['usgs_named']+r['usgs_generic']:
            a=read(e)
            if a:
                a=dict(a,match=e['match'],condition=a.get('condition') or ('weathered_or_coated' if re.search(r'weather|alt\.|varnish|jarosite|hematitic|hydrated|goethite',e['title'],re.I) else 'unspecified'))
                rows.append(a)
        unique={x['sample_id']:x for x in rows};rows=list(unique.values())
        rocks[r['id']]=dict(observations=rows,statistics=stats(rows),qualification='USGS_named_or_generic_reference_not_a_universal_rock_population' if rows else 'no_usable_USGS_rock_reference_keep_mineral_model_no_fitted_rock_tint')
    # Physical mixtures verified against original USGS specimen descriptions.
    recipes=[(12536,{'hematite':.02,'quartz':.98},'hema_qtz_gds76',125e-6),
             (12463,{'goethite':.02,'quartz':.98},'goethite_gds240',50e-6),
             (12075,{'calcite':.8,'kaolinite':.2},'calcite_mixture_gds213',50e-6)]
    fits=[]
    for rec,mode,doc,d in recipes:
        a=read(next(x for x in idx if x['record']==rec));assert a, (rec, rejects[-1])
        fits.append(dict(kind='intimate_mixture',observation=a,mode=mode,grain_m=d,grain_qualification='nuisance prior: ground sample size unknown; GDS76 upper250um midpoint; others50um; sensitivity required',metadata='https://pubs.usgs.gov/of/2003/ofr-03-395/DESCRIPT/S/'+doc+'.html'))
    for rec,realization in [(11848,'basalt_tholeiitic.reference'),(2572,'chert.white')]:
        a=read(next(x for x in idx if x['record']==rec));assert a, (rec, rejects[-1])
        fits.append(dict(kind='whole_rock_proxy',observation=a,realization=realization,qualification='Measured rock; catalogue composition/size NOT specimen-matched; weak proxy fit weight0.25. No slab/particulate inferred.'))
    # Coating targets with unknown thickness are not identifiable, so retain v4 shell priors.
    obj=dict(schema='pf-rockspec-optics/5',composition_schema='pf-rockspec/3',method=METHOD,source='USGS only; Kokaly et al. 2017 doi:10.5066/F7RR1WDJ',minerals=minerals,rocks=rocks,transport_transition_lengths=[1.,5.],oxide_scattering_multiplier=1.,transport_fit_qualification='pending_USGS_fit',coating_qualification='v4 stock-conserving shell priors retained: no specimen-matched substrate, coverage and thickness calibration',missing_minerals=missing)
    if len(sys.argv)>5:
        fitpath=Path(sys.argv[5]);fit=json.load(open(fitpath));best=fit['best']
        obj['transport_transition_lengths']=best['limits']
        obj['oxide_scattering_multiplier']=best['oxide_scattering_multiplier']
        obj['transport_fit_qualification']='USGS-only 5 targets:3 physical mixtures with documented endmembers;2 weak whole-rock modal proxies. 324-point conditional grid. Weathered limestone excluded. Grain sizes and optical mapping remain priors; boundary optimum not universal identification.'
        obj['fit_receipt']=dict(best=best,grid_count=fit['grid_count'],sha256=hashlib.sha256(fitpath.read_bytes()).hexdigest())
    out.write_text(json.dumps(obj,indent=2)+'\n')
    (lane/'fit_targets_usgs.json').write_text(json.dumps(fits,indent=2)+'\n')
    (lane/'extraction_audit.json').write_text(json.dumps(dict(method=METHOD,missing=missing,accepted_mineral_observations=sum(len(m['observations']) for m in minerals.values()),rejections=rejects,grain_known=sum(x['grain_m'] is not None for m in minerals.values() for x in m['observations']),input_sha256={str(covdir/f):hashlib.sha256((covdir/f).read_bytes()).hexdigest() for f in ['coverage.json','splib07_index.json']},wavelength_sha256=hashlib.sha256(wlfile.read_bytes()).hexdigest()),indent=2)+'\n')
    print('minerals',sum(bool(x['observations']) for x in minerals.values()),'observations',sum(len(x['observations']) for x in minerals.values()),'rocks',sum(bool(x['observations']) for x in rocks.values()),'missing',missing)
if __name__=='__main__':main()
