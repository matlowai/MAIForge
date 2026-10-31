#!/usr/bin/env python3
"""Offline morphology priors. No legacy catalogue writes; no fitted rock colours."""
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
c = json.loads((root / 'data/catalog_v3.json').read_text())

def prior(v, lo, hi, unit, conditions):
    return dict(value=v, interval=[lo, hi], unit=unit,
                interval_kind='plausible_prior_support', distribution='unspecified_joint_prior',
                reference='OPT4', evidence_kind='project_prior',
                support_scale='conditional sediment morphology; not measured',
                conditions=conditions, covariance_group='mass_constrained_shell')

profiles = {}
for r in c['realizations']:
    if r['lithofacies'] not in ['sandstone_quartz', 'arkose', 'greywacke', 'conglomerate', 'sand', 'regolith']:
        continue
    modes = {p['mineral']: p['mass_fraction']['value'] for p in r['spec']['mode']}
    # Partition existing iron stock, never add pigment or invent hematite.
    iron = modes['goethite']
    cemented = r['lithofacies'] in ['sandstone_quartz', 'arkose', 'greywacke', 'conglomerate']
    coating = iron * .8 if cemented else iron
    cement = ([dict(mineral='quartz', mass_fraction=.10),
               dict(mineral='goethite', mass_fraction=iron-coating)] if cemented else [])
    profiles[r['id']] = dict(
        coating_mode=[dict(mineral='goethite', mass_fraction=coating)],
        cement_mode=cement,
        coverage=prior(.7, .2, .95, '1', 'patchy diagenetic grain coats; prior shared by all sediments'),
        shell_thickness_support_m=prior(1e-6, 1e-8, 5e-5, 'm', 'nominal scale only; actual thickness solved from stock, coverage and physical grain diameter'),
        pigment_scattering_length_m=prior(1e-6, 1e-7, 1e-5, 'm', 'fine Fe oxyhydroxide aggregate effective scattering length; not measured optical constants'),
        cement_grain_m=prior(5e-6, 1e-6, 5e-5, 'm', 'microcrystalline silica/iron cement; no cement in loose sand or regolith'),
        qualification='Existing dry mineral mass partitioned into grain cores, cement and coats. Lithic source identities remain unresolved. Goethite used because admitted modes contain it; no invented hematite.')

out = dict(schema='pf-rockspec-optics/4', composition_schema='pf-rockspec/3',
           references={'OPT4': dict(
               citation='Kubelka and Munk (1931), two-flux layer theory; Hapke (2012), Theory of Reflectance and Emittance Spectroscopy, 2nd ed.; Cornell and Schwertmann (2003), The Iron Oxides, 2nd ed.; Tucker (2001), Sedimentary Petrology, 3rd ed.',
               qualification='Offline bibliographic recollection only. References support mechanisms, not the numerical morphology priors. No pages, quotes, spectral extraction or library calibration claimed.',
               verification='offline_unverified_recollection')},
           transport_transition_lengths=[1., 5.],
           transition_qualification='Reduced two-flux diffuse penetration depth. Smoothstep from one to five return-weighted penetration lengths; global model choice, not a measured universal transition.',
           sediments=profiles)
(root / 'data/optics_v4.json').write_text(json.dumps(out, indent=2)+'\n')
print(f'{len(profiles)} sediment morphology profiles; no v1/v2/v3 file changed')
