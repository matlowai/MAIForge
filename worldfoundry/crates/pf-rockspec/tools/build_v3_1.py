#!/usr/bin/env python3
"""Offline v3.1b data/fixture generator. Reads frozen v3/v2, writes only v3_1.
References are recalled books, not consulted extracts or measured specimen data.
"""
import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
c = json.loads((ROOT / 'data/catalog_v3.json').read_text())
c['schema_version'] = 'pf-rockspec-assemblage/3.1b'
c['core']['schema_version'] = 'pf-rockspec/3.1b'
c['core']['numerical_model'] = 'assemblage+connected-km/3.1b'
c['core']['references']['RS31'] = dict(
    citation='Winter (2010), Principles of Igneous and Metamorphic Petrology, 2nd ed.; '
             'Deer, Howie and Zussman (1992), An Introduction to the Rock-Forming Minerals, 2nd ed.; '
             'Butt and Zeegers, eds. (1992), Regolith Exploration Geochemistry in Tropical and Subtropical Terrains; '
             'MODEL_V3_1.md',
    qualification='Offline memory only; broad petrological relationships, no consulted pages or extracted modal/optical numbers. '
                  'All new numbers are conditional project priors, not setting-wide averages, measurements, or calibration.',
    verification='offline_memory_unverified_project_prior')

def q(v, lo, hi, unit='1', conditions='conditional scenario; not a setting-wide distribution'):
    return dict(value=v, interval=[lo, hi], unit=unit,
                interval_kind='plausible_prior_support', distribution='unspecified_joint_prior',
                reference='RS31', evidence_kind='project_prior',
                support_scale='body realization; no site admission', conditions=conditions,
                covariance_group='v31_closed_assemblage')

def mode(parts, conditions='dry mass mode; joint simplex; illustrative differentiated or olivine-bearing lava interior'):
    assert abs(sum(parts.values())-1) < 1e-12
    return [dict(mineral=k, mass_fraction=q(v, max(0, v*.3), min(1, v*1.8),
             conditions=conditions))
            for k, v in parts.items()]

minerals = {m['id']: m for m in c['core']['minerals']}
real = {r['id']: r for r in c['realizations']}
realization_priors = {}
# Specific permissible members of overlapping basalt populations, never a claim
# that tectonic setting uniquely determines mode. Reference modes stay unchanged.
for rid, parts, cond in [
    ('basalt_tholeiitic.lip', dict(plagioclase=.50, augite=.34, orthopyroxene=.02,
        olivine=.02, magnetite=.06, volcanic_glass=.06),
        'declared moderately differentiated flood-lava assemblage; not every LIP basalt'),
    ('basalt_tholeiitic.island_shield', dict(plagioclase=.45, augite=.32,
        olivine=.12, magnetite=.04, volcanic_glass=.07),
        'declared olivine-bearing shield-lava assemblage; not every ocean-island tholeiite'),
]:
    r = real[rid]; s = r['spec']
    s['mode'] = mode(parts)
    realization_priors[rid] = dict(
        scenario=cond, mode=copy.deepcopy(s['mode']),
        qualification='Modal change only. All grain, porosity, An and structure parameters retain v3 values. Offline illustrative prior, not extracted measurement.')
    rho = (1-s['intrinsic_porosity']['value']) / sum(
        v / minerals[k]['density_kg_m3']['value'] for k, v in parts.items())
    s['dry_bulk_density_kg_m3'] = q(rho, rho*.8, rho*1.2, 'kg/m3', 'derived nominal density; resolver applies body An')
    dep = rho * (1-s['interparticle_porosity']['value'])
    s['dry_deposit_density_kg_m3'] = q(dep, dep*.8, dep*1.2, 'kg/m3', 'derived nominal deposit density')
    s['crystallinity'] = q(1-parts['volcanic_glass'], 0, 1)
# Alkali reference and cap retain their complete v3 records.

# Explicit wet-temperate residual surface scenarios. No kinetics or bulk mutation.
# BIF has no Al silicates in this mode; generic silicate hydrolysis was erroneous.
residues = {
    'pegmatite': dict(kaolinite=.8, quartz=.2),
    'skarn': dict(goethite=.45, smectite=.35, calcite=.2),
    'banded_iron_formation': dict(goethite=.5, hematite=.25, quartz=.25),
    'laterite': dict(goethite=.4, hematite=.15, gibbsite=.3, kaolinite=.15),
    'mylonite': dict(kaolinite=.65, quartz=.25, goethite=.1),
    'gossan': dict(goethite=.65, hematite=.15, quartz=.2),
}
for rock in c['core']['rocks'] + [r['spec'] for r in c['realizations']]:
    if rock['id'] in residues:
        rock['weathering_products'] = list(residues[rock['id']])
        if rock['id'] == 'banded_iron_formation':
            rock['weathering_family'] = 'oxidation_residual'

old_surfaces = {}
for row in json.loads((ROOT/'data/swatch_inputs_v2.json').read_text())['samples']:
    old_surfaces[(row['material_ref']['spec_id'], row['surface']['surface_id'].split('/')[-1])] = row['surface']
surfaces = []
weather_evidence = {}
for r in c['realizations']:
    key = r['lithofacies']; family = r['spec']['weathering_family']
    if key not in residues:
        continue
    for history in ['fresh_dry', 'fresh_wet', 'weathered_dry']:
        s = copy.deepcopy(old_surfaces[('granite', 'fresh_dry')])
        s['surface_id'] = key+'/'+history
        if history == 'fresh_wet':
            s.update(wet_coverage=1., water_film_m=.0001,
                     forcing_record='manufactured fully wet fresh fracture; 100 micrometre film')
        if history == 'weathered_dry':
            coverage = .18 if key in ['pegmatite','mylonite'] else .4
            thickness = .00005
            s.update(mechanism='inherited', exposure_elapsed_s=None,
                creation_event='diagnostic-inherited-v31', rms_slope_history_rad=.08, abrasion=.05,
                forcing_record=f'manufactured wet-temperate {family} residue on {key}; supplied surface stock; no kinetics, elapsed age, or bulk reaction claimed',
                coating=dict(source_event=f'manufactured-v31-{key}-residue',
                    mineral_mode=[dict(mineral=k,mass_fraction=v) for k,v in residues[key].items()],
                    grain_diameter_m=.0002, thickness_m=thickness, coverage=coverage,
                    extinction_m_inv=c['core']['model']['rind_extinction_m_inv']['value']))
            weather_evidence[key] = dict(family=family,
                residue_mode=mode(residues[key], 'supplied diagnostic residue dry mass mode; joint simplex; no chemical reaction or external flux solved'),
                coverage=q(coverage, .05, .65, conditions='manufactured wet-temperate residual patch coverage, not predicted exposure'),
                thickness_m=q(thickness, .00001, .0002, 'm', 'manufactured residual rind; unknown weathering duration'),
                grain_diameter_m=q(.0002,.00002,.001,'m','manufactured residual aggregate domain'),
                rms_slope_history_rad=q(.08,.02,.2,'rad','manufactured surface roughening'),
                abrasion=q(.05,0,.2,conditions='manufactured minor wear'),
                qualification='Supplied diagnostic surface phase assemblage. No reaction, elemental balance, or stock transfer inferred.')
        for footprint in [.1, 2., 50.]:
            surfaces.append(dict(id=f'{r["id"]}/{history}/{footprint}', realization_id=r['id'],
                history=history, weathering_family=family, surface=s, footprint_m=footprint))

def write(name, data):
    (ROOT/'data'/name).write_text(json.dumps(data, indent=2)+'\n')
write('catalog_v3_1.json', c)
write('realization_priors_v3_1.json', dict(reference=c['core']['references']['RS31'], realizations=realization_priors))
qual = 'Manufactured conditional v3.1b surfaces for six new families only; evaluated by unchanged v3 optics; no world admission or reaction rates'
write('swatch_inputs_v3_1.json', dict(schema='pf-rockspec-surfaces/3.1b', qualification=qual, samples=surfaces))
write('diagnostic_inputs_v3_1.json', dict(schema='pf-rockspec-surfaces/3.1b', qualification=qual,
    samples=[s for s in surfaces if s['realization_id'].split('.')[0] in residues]))
write('weathering_priors_v3_1.json', dict(reference=c['core']['references']['RS31'], families=weather_evidence))
print(f'v3.1b: {len(real)} realizations; {len(surfaces)} swatches; 54 new-family diagnostics')
