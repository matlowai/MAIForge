# pf-rockspec v3.1b: frozen weathering surfaces, v3 optics

This replaces the unmerged v3.1 candidate after the integrator's 2026-10-25
04:37 EDT rejection. The areal/transport switch made coarse granite and basalt
unrealistically bright. Per-grain mineral reflectance/scattering priors remain
uncalibrated, so correcting mixing alone did not deliver physical realism.
The rejected source is preserved as a patch in the lane's `preserved-v31/`;
its original reports, images and exports remain evidence, not active guidance.

File and Rust module names retain `v3_1` for this replacement. The data identity
is explicitly **pf-rockspec/3.1b**, catalogue schema
**pf-rockspec-assemblage/3.1b**, and model **assemblage+connected-km/3.1b**.
Use `pf_rockspec::v3_1::catalogue()?`, then supply the ordinary formation and body
admission with a /3.1b MaterialRef. Old default constructors remain v1/v2/v3.
Rejected /3.1 data/model combinations are no longer accepted. V4 continues to
accept only its existing /3 composition contract.

No optical equation, mineral reflectance, scattering prior, roughness or Fresnel
formula is changed. The complete v3 optical path is shared by v3.1b. There is
no grain-transport/areal switch in this revision and no v4 optical call in its
exporter. V1/v2/v3/v4 data and exporters remain unchanged.

This crate is the engine's shared body/fracture/render material subsystem. The
new surfaces are explicit catalogue diagnostics. This arc does not change a
world producer, production surface history, inventory admission or Unreal consumer.

## Granite: known limitation, not fixed

K-feldspar colour already reaches the v3 optical path. Both I/S pairs contain
0.32 K-feldspar by dry mass, about 0.3322/0.3300 by solid volume. White resolves
to aggregate linear RGB [0.74, 0.73, 0.70]; pink retains [0.64, 0.47, 0.39].
In the inherited intimate Kubelka-Munk mixture, the 10% biotite/hornblende mass
fraction supplies 91.2 to 96.1% of RGB absorption/scattering ratio K/S for I-type
pink. The 4% biotite supplies 82.2 to 88.8% for S-type pink. These are computed
model contributions, not measured photon paths. That nonlinear mixture suppresses
the feldspar contrast. The colour parameter is not disconnected or lost.

The rejected areal model exposed a second problem: uncalibrated per-grain priors
made coarse rocks too bright. V3.1b therefore leaves all granite data and optics
unchanged. The maximum pink/white linear-channel differences at 2 m remain
2.1083/255 for I-type and 5.3269/255 for S-type. Tests require all 54 granite
surface/scale outputs to match v3 exactly, apart from version metadata. No
brightness, saturation or tint correction is applied. Calibration remains pending
model specification on a spectral library; no acquisition or calibration job is launched.

## Conditional basalt modal differences only

The original v3 tholeiitic LIP and shield modes were identical. Only these two
realizations receive new closed dry-mass modes:

| Realization | Modal dry mass fractions |
|---|---|
| basalt_tholeiitic.lip | plagioclase .50, augite .34, orthopyroxene .02, olivine .02, magnetite .06, volcanic_glass .06 |
| basalt_tholeiitic.island_shield | plagioclase .45, augite .32, olivine .12, magnetite .04, volcanic_glass .07 |

These are illustrative differentiated flood-lava and olivine-bearing shield-lava
assemblages. They are permissible members of overlapping populations, not a
claim that tectonic setting determines a unique mode. All grain sizes, An values,
porosity, structure, mechanics and surface parameters retain v3 values. Density
and crystallinity are recomputed solely from the changed mode. The reference
tholeiite and both alkali basalt realizations are completely unchanged.

The modes are retained because the measured fresh linear luminance Y remains
within the requested 0.05 to 0.10 interval at 0.1, 2 and 50 m. Across all five
basalt realizations and three footprints, Y spans 0.05087 to 0.05865. At 2 m,
LIP changes 0.05114 to 0.05222 and shield 0.05114 to 0.05418. These subtle optical
differences do not establish visually separable tectonic settings. The alkali
reference/cap pair remains nearly indistinguishable under v3 at these footprints.
No values were tuned to a brightness target. The already proposed modes pass the
requested range with the unchanged v3 closure.

`realization_priors_v3_1.json` records the two conditional mode priors. Each
fraction has a broad support, unit, reference and conditions; fractions live on
a joint simplex and must not be independently sampled. Offline recalled
petrology supports possible assemblages, not these exact measured percentages.

## Six frozen surface families

Both `swatch_inputs_v3_1.json` and `diagnostic_inputs_v3_1.json` now contain
exactly **54 cases**: six families x fresh dry/fresh wet/weathered dry x
0.1/2/50 m. Their surface records are retained from the rejected candidate;
they are evaluated with the unchanged v3 optical path. The exporter additionally
assembles a 603-case comparison grid using existing v2 surface records for other
families. Those inherited surfaces are not frozen again in new data files.

| Family | Weathering family | Supplied residue dry mass mode |
|---|---|---|
| pegmatite | silicate_hydrolysis | kaolinite .80, quartz .20 |
| skarn | silicate_hydrolysis | goethite .45, smectite .35, calcite .20 |
| banded_iron_formation | oxidation_residual | goethite .50, hematite .25, quartz .25 |
| laterite | oxidation_residual | goethite .40, hematite .15, gibbsite .30, kaolinite .15 |
| mylonite | silicate_hydrolysis | kaolinite .65, quartz .25, goethite .10 |
| gossan | oxidation_residual | goethite .65, hematite .15, quartz .20 |

BIF's Al-free quartz/iron-oxide mode does not support the old generic silicate
hydrolysis label. Only the new catalogue corrects it to oxidation/residual.
Weathering-product lists in v3.1b describe these supplied residues, including
resistant surviving phases. BIF receives no invented Al clay and pegmatite receives
no invented iron stain. Skarn can retain calcite alongside products of Fe-silicate
alteration, but no dissolution/oxidation rate is inferred.

The explicit wet-temperate surfaces have 50 micrometre rind thickness,
200 micrometre aggregate domain, 0.08 rad additional slope and .05 abrasion.
Coverage is .18 for pegmatite/mylonite and .40 for the other four families.
`weathering_priors_v3_1.json` supplies supports, units and conditions. Fresh wet
has a 100 micrometre water film and no rind. Weathered exposure elapsed time is
null. The evaluator never creates a surface or alters bulk minerals. Residues are
manufactured supplied surface stocks, not solved reactions or world-stock transfers.

Tests compare all 54 cases against direct v3 evaluation of the same surface.
They validate product families, residue mass closure, exact prior/fixture agreement,
absence of implicit fresh-face rinds, passivity and unchanged body state. No
historical v3 weather fixture existed for these six families: reports label that
baseline absent. A same-surface v3 evaluation is a control, not historical evidence.

## Offline evidence and validation

RS31 cites Winter (2010), *Principles of Igneous and Metamorphic Petrology*,
second edition; Deer, Howie and Zussman (1992), *An Introduction to the
Rock-Forming Minerals*, second edition; and Butt and Zeegers, editors (1992),
*Regolith Exploration Geochemistry in Tropical and Subtropical Terrains*.
All are from offline memory only. No pages were consulted. New numerical modes
and residues are project priors, not extracted measurements or confidence intervals.
No network or spectral library was used.

Run `cargo test --offline --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml`
with the assigned `CARGO_TARGET_DIR` and `CARGO_BUILD_JOBS=8`. Generate data with
`tools/build_v3_1.py`; export using `rockspec_export_v3_1 OUTPUT_DIR`; render tables
and sheets with `tools/report_v3_1.py OUTPUT_DIR REPORT_DIR` using existing Pillow.
Tests do not read lane files. All old default exports, including images, are rerun
and SHA256-compared to the original pre-lane baseline.

This revision does not claim calibrated granite colour, basalt setting recognition,
weathering kinetics, chemical conservation, world admission or game-visible
integration. The next work is integrator review of these bounded fixtures and
modal changes. Optical calibration remains a separate pending decision.
