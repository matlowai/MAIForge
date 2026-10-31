#!/usr/bin/env python3
"""Audit actual exports and retained baseline. No material mutation or data fitting."""
import csv, hashlib, json, re, sys
from pathlib import Path

fig, ev, heavy = map(Path, sys.argv[1:4])
root = Path(__file__).resolve().parents[1]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
def write(p, obj): p.write_text(json.dumps(obj, indent=2)+'\n')
lum = lambda out: sum(a*b for a,b in zip(out['base_reflectance_linear'], [.2126,.7152,.0722]))
rows = json.loads((fig/'samples_v4.json').read_text())['samples']
diag = json.loads((fig/'diagnostics_v4.json').read_text())['samples']
assert len(rows)==67 and len(diag)==128

baseline=json.loads((ev/'baseline_files.json').read_text())
changed=[p for p,h in baseline.items() if sha(Path(p))!=h]
assert changed==['worldfoundry/crates/pf-rockspec/src/lib.rs','worldfoundry/crates/pf-rockspec/src/v3.rs'],changed
data=[dict(path=p,sha256=h) for p,h in baseline.items() if '/data/' in p]
assert all(sha(Path(x['path']))==x['sha256'] for x in data)
exports=[]
for p in sorted((heavy/'before').iterdir()):
    if not p.is_file(): continue
    q=heavy/'after'/p.name
    assert p.read_bytes()==q.read_bytes(),p.name
    exports.append(dict(file=p.name,bytes=p.stat().st_size,sha256=sha(p)))
write(ev/'keyoff_identity.json',dict(legacy_data=data,unchanged_original_files=len(baseline)-len(changed),
    intentional_source_hooks=changed,byte_identical_exports=exports))

# These are review aids, never tests/calibration targets or fitted model inputs.
bands={
 'basalt_tholeiitic':(.05,.10,'existing operator sanity prior'),
 'granite':(.25,.35,'existing operator sanity prior'),
 'marble':(.60,.90,'existing operator sanity prior'),
 'obsidian':(.008,.03,'existing preservation envelope'),
 'andesite':(.12,.25,'design review F3/F6 offline recollection'),
 'dacite':(.25,.40,'design review F3/F6 offline recollection'),
 'diorite':(.20,.30,'design review F3/F6 offline recollection'),
 'schist':(.20,.35,'design review F6 muscovite-schist recollection'),
 'peridotite':(.15,.30,'design review F6 fresh olive-green recollection'),
 'scoria':(.04,.10,'design review F6 recollection'),
 'sandstone_quartz':(.25,.40,'design review F6 buff sandstone recollection'),
 'siltstone':(.20,.35,'design review F6 recollection'),
 'chalk':(.80,.90,'design review F6 recollection'),
 'coal':(.03,.06,'design review F6 recollection'),
 'pumice':(.50,.70,'design review F6 recollection'),
 'rhyolite':(.30,.50,'design review F6 recollection'),
 'limestone':(.30,.60,'design review F6 recollection'),
 'silt_clay':(.20,.35,'design review F6 recollection'),
}
special={
 'shale.organic_rich':(.04,.08,'design review F6 black shale recollection'),
 'shale.grey':(.10,.25,'design review F6 grey shale recollection'),
 'chert.flint':(.05,.15,'design review F6 dark flint recollection'),
 'chert.white':(.50,.70,'design review F6 white chert recollection'),
 'ignimbrite.unwelded':(.50,.70,'design review F6 pale unwelded ash recollection'),
 'ignimbrite.welded':(.15,.35,'design review F6 welded ignimbrite recollection'),
}
comments={
 'granodiorite':'Coarse dark bias removed, but Y above 0.6 is visibly too pale for a representative mixed hand specimen. No calibrated numerical band.',
 'gabbro':'Coarse mix is visibly much too pale for the intended dark mafic specimen; mineral scattering priors remain unresolved.',
 'amphibolite':'Hornblende/plagioclase contrast becomes visible; diffuse brightness likely high. No adopted numerical band.',
 'gneiss':'Band contrast is visible, but felsic mean is very bright. No calibrated phase/size populations.',
 'arkose':'Iron coating mechanism darkens the pale v3 core. Buff/pink appearance is plausible, not calibrated.',
 'greywacke':'Matrix/iron-bearing mixture stays darker than coated quartz sandstone. No qualified specimen band.',
 'conglomerate':'Coated mineral-domain proxy is plausible; actual lithic populations and clast topology remain unresolved.',
 'sand':'Coated sand is darker than the same iron stock dispersed. No universal sand-albedo target.',
 'regolith':'Coated residual mixture darkens relative to dispersed stock. Actual weathering history and site colour unknown.',
 'basalt_alkali':'Brightening is required by the supplied grain scale in this closure. The island-cap realization is too pale for a typical basalt; not accepted.',
 'ignimbrite':'Reference is explicitly pale unwelded ash. Still bright; do not treat it as a welded-ignimbrite success.',
 'peridotite':'Fresh green channel exceeds red and blue; near the high edge of the unverified fresh range.',
 'obsidian':'Dark glass preserved; tiny modal closure change is exposed.',
 'marble':'Near-monomineralic response preserved; tiny areal/intimate difference is exposed.',
 'chalk':'Fine-grained intimate limit preserved.',
 'shale':'Fine-grained intimate limit preserved; organic-rich variant stays dark.',
}
verdicts=[]; table=[]; max_energy=0.; min_reflectance=1.; max_reflectance=0.
for r in rows:
    s=r['scales'][1]; value=lum(s['v4']); band=special.get(r['id'],bands.get(r['facies']))
    status='unqualified: no justified offline brightness interval'
    if band:
        status=('too dark' if value<band[0] else 'too bright' if value>band[1] else 'within unverified sanity interval')
    verdicts.append(dict(id=r['id'],v3_Y=lum(s['v3']),v4_Y=value,status=status,
        interval=None if not band else band[:2],basis=None if not band else band[2],
        note=comments.get(r['facies'],'A conditional material appearance only; no specimen or spectral validation claimed.')))
    for s in r['scales']:
        table.append(dict(id=r['id'],footprint_m=s['footprint_m'],v3_Y=lum(s['v3']),v4_Y=lum(s['v4']),
            delta_Y=lum(s['v4'])-lum(s['v3']),v4_dispersed_coat_Y=lum(s['v4_without_coatings'])))
        for v in ['v3','v4','v4_without_coatings']:
            out=s[v]
            for j in range(3):
                max_energy=max(max_energy,out['diffuse_energy_linear'][j]+out['f0_linear'][j])
                min_reflectance=min(min_reflectance,out['base_reflectance_linear'][j]);max_reflectance=max(max_reflectance,out['base_reflectance_linear'][j])
                assert 0<=out['base_reflectance_linear'][j]<=1
                assert out['diffuse_energy_linear'][j]+out['f0_linear'][j]<=1+1e-12
for r in diag:
    for j in range(3):
        out=r['output'];assert 0<=out['base_reflectance_linear'][j]<=1
        assert out['diffuse_energy_linear'][j]+out['f0_linear'][j]<=1+1e-12
        max_energy=max(max_energy,out['diffuse_energy_linear'][j]+out['f0_linear'][j])
write(ev/'optical_comparison_v4.json',table)
write(ev/'per_realization_verdicts.json',verdicts)
with (ev/'optical_comparison_v4.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=table[0].keys());w.writeheader();w.writerows(table)

lines=['# Every realization: v3 versus v4 linear luminance Y','',
       'Same body composition and surface. Three footprint scales in metres. No fitted target or spectral calibration. Full precision and coating ablation in evidence/optical_comparison_v4.csv.', '',
       '| Realization | v3 0.1m | v4 0.1m | v3 2m | v4 2m | v3 50m | v4 50m |',
       '|---|---:|---:|---:|---:|---:|---:|']
for r in rows:
    nums=[lum(s[v]) for s in r['scales'] for v in ['v3','v4']]
    lines.append('| '+r['id']+' | '+' | '.join(f'{n:.6f}' for n in nums)+' |')
lines+=['','## Coat morphology ablation at 2 m','',
        'Without coating means the same coat stock is dispersed into the cores. No mass is deleted. Cement stays partitioned.', '',
        '| Realization | v3 Y | v4 dispersed Y | v4 coated Y | Thickness um | Coverage | Cement dry mass |',
        '|---|---:|---:|---:|---:|---:|---:|']
for r in rows:
    m=r['morphology']
    if m:
        s=r['scales'][1];lines.append(f"| {r['id']} | {lum(s['v3']):.6f} | {lum(s['v4_without_coatings']):.6f} | {lum(s['v4']):.6f} | {m['thickness_m']*1e6:.6f} | {m['coverage']:.2f} | {m['cement_mass_fraction']:.6f} |")
lines+=['','## Honest per-realization disposition','',
        'Intervals are inherited operator sanity checks or design review offline recollections, not measurements. A value inside a band does not certify appearance. Coarse-rock overshoots are failures for default visual adoption, not hidden as model passes.', '',
        '| Realization | 2 m Y | Verdict | Basis and caveat |','|---|---:|---|---|']
for v in verdicts:
    basis=('No calibrated interval.' if v['interval'] is None else f"{v['interval']}; {v['basis']}.")
    lines.append(f"| {v['id']} | {v['v4_Y']:.6f} | {v['status']} | {basis} {v['note']} |")
(ev.parent/'CHANGED_NUMBERS.md').write_text('\n'.join(lines)+'\n')
tests=(ev/'tests_third.log').read_text();counts=[int(n) for n in re.findall(r'test result: ok\. (\d+) passed',tests)]
assert sum(counts)==48
write(ev/'acceptance.json',dict(realizations=67,footprints=[.1,2,50],comparison_rows=len(table),
    diagnostics=len(diag),sediment_profiles=sum(r['morphology'] is not None for r in rows),
    tests_passed=sum(counts),test_groups=counts,legacy_data_files=len(data),identical_export_files=len(exports),
    reflectance_min=min_reflectance,reflectance_max=max_reflectance,max_diffuse_plus_f0=max_energy,
    calibration_gate='NOT CLAIMED',native_render_gate='NOT CLAIMED',
    coarse_default_visual_acceptance='FAIL: several coarse rocks too bright; retain opt-in pending calibration'))
print(json.dumps(json.loads((ev/'acceptance.json').read_text()),indent=2))
