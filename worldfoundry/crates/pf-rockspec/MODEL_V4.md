# pf-rockspec optics v4

Opt in with `v4::Optics::embedded()?.surface_sample(&resolved_v3_body, &surface, footprint)`.
Default v1, v2 and v3 evaluators are unchanged. The body, its `MaterialRef`, admitted
mineral stock, formation, mechanics and density stay v3. The output optical contract
is `pf-rockspec-optics/4`, model `grain-transport+finite-km-shell/4`. This is an
explicit optical evaluation version, not a new composition version or a silent
migration. `data/optics_v4.json` contains the additional morphology priors.

This belongs to the engine's shared geological material subsystem, consumed by
body/fracture/render adapters. It uses the actual shared CPU surface evaluator and
v3 resolved bodies. The lane cannot modify the world/inventory/Unreal consumers.
The included images and 128 cases are manufactured material diagnostics, not a
new standalone application or evidence of in-game adoption. Native integration
and measured spectral calibration are not claimed.

## Grain transport

For each RGB channel and mineral i, retain the inherited aggregate reflectance
R0, reference diameter d0 and internal scattering length l. Infer effective
absorption and scattering coefficients, in inverse metres:

```
q0_i = (1-R0_i)^2 / (2 R0_i)
K_i  = q0_i (1/d0_i + 1/l_i)
S_i  = 1/d_eff + 1/l_i
d_eff = d / optical_boundary_efficiency
R_i  = 1 / (1+K_i/S_i + sqrt((K_i/S_i)^2 + 2 K_i/S_i))
v_i  = (w_i/rho_i) / sum_j(w_j/rho_j)
```

The random-section area fraction is volume fraction v_i, not dry mass fraction
w_i. The physical mineral optical domain d is `optical_domain_diameter_m`, before
boundary-efficiency adjustment. A conglomerate clast can contain many such mineral
grains; centimetre clast d50 must not become the mineral transport diameter.
Each current realization has one monodisperse optical domain approximation. A
per-species size distribution and explicit lithic-fragment topology remain absent.

Two limiting responses are evaluated:

```
R_int = KM(sum(v_i K_i) / sum(v_i S_i))
R_area = sum(v_i R_i)
ell_i = 1/sqrt(K_i (K_i+2 S_i))
ell_return = sum(v_i R_i ell_i) / sum(v_i R_i)
x = clamp((d/ell_return - 1)/4, 0, 1)
a = x*x*(3-2*x)
R = (1-a)*R_int + a*R_area
```

ell_i is the diffuse penetration/attenuation scale in the two-flux equations,
not a measured crystal mean free path. Return weighting makes the bright grains
that carry surviving photons govern the lateral sampling approximation; rare
opaque grains cannot force an otherwise transparent mixture into the areal limit.
The physical hypothesis is that a photon returning within one grain samples fewer
other constituents when d spans several penetration lengths. One to five lengths
is a global reduced-order transition with zero endpoint slopes. It is a stated
model choice, not a fitted law or a rock-key brightness table. RGB transitions can
differ because absorption differs. Footprint, camera, rock key and target Y never
enter the switch. Boundary efficiency affects S and hence transport, while d
remains the actual mineral-domain scale.

For exact K=0, penetration length is infinite, encoded as f64::MAX in the numeric
receipt; the transition remains intimate. Current catalogue priors all absorb.
All mixture weights are nonnegative and sum to one. Each endpoint is in [0,1], so
the blend cannot create energy. For fixed equal S, KM convexity makes the areal
endpoint at least as bright as the intimate endpoint. Increasing the closure
weight then monotonically brightens a heterogeneous mixture. Increasing d also
reduces boundary scattering and darkens an individual mineral. These competing
effects mean that a mixed rock has no universal monotone Y versus grain-size law.
Tests separate these cases instead of imposing an arbitrary monotonicity checker.

**An additional consequential correction:** v3 averaged K_i/S_i. V4 averages K and
S separately for the intimate limit. They agree when S is common, which covers
most inherited minerals. They differ for welded glass with its longer internal
scattering length. The welded ignimbrite change is therefore real and reported;
it is not hidden behind a rock-specific exception.

## Grain coatings and cement

Six explicit conditional profiles cover quartz sandstone, arkose, greywacke,
conglomerate, sand and regolith. They partition existing dry mass into core,
cement and coating. No mineral is added, no biotite hack restored, no density
changed and no ore admitted. The closed mode remains the authority.

For cemented profiles, 80% of the existing goethite stock coats grains and 20%
is oxide cement; 0.10 of dry mass is silica cement taken from existing quartz.
For sand and regolith, all existing goethite is grain coating, with zero cement.
These are new **conditional morphology priors**, not discoveries about a seeded
world. The v3 sediment parameters remain historical v3 values. V4 deliberately
repartitions their same stock, recorded explicitly in `optics_v4.json`.

Goethite is used because these admitted modes contain goethite. Hematite is not
invented to achieve a colour. The equations accept a supplied hematite partition
if a body's admitted mode funds it. Existing lithic-derived mineral constituents
are retained, but no source-rock/lithic geometry is inferred from mineral modes.

Coverage c = 0.70 has support [0.20,0.95]. For a spherical core of diameter d, the
coat thickness t is solved from core and coating volumes, not selected to hit Y:

```
V_f = c V_core ((1+2t/d)^3 - 1)
t = d/2 ((1 + V_f/(c V_core))^(1/3) - 1)
```

Volumes are m3 per kg dry solid. Cement is excluded from V_core. This closes the
shell volume exactly; partial coverage is a patchy shell approximation. A 0.01 to
50 micrometre thickness support is recorded, with 1 micrometre as the nominal
prior scale, not a value used to override the stock-constrained solution. The
export reports the actual thickness and whether it lies in that support.

The pigment aggregate has an effective scattering-length prior of 1 micrometre
with support [0.1,10] micrometres. Its K/S is inferred from the existing pigment
RGB R0, while S = 1/l_pigment sets the dimensional scale. Thus this is a new,
explicit fine-aggregate coefficient prior, **not measured Fe-oxide optical
constants**. Cement has a 5 micrometre domain, support [1,50] micrometres, and
uses the intimate mixture of its constituents. Coated cores and cement responses
are mixed by their solid volume/random-section fractions.

For a homogeneous reciprocal two-flux slab, gamma=sqrt(K(K+2S)), e=exp(-gamma t),
u=1-exp(-2 gamma t), D=gamma(1+e^2)+(K+S)u:

```
R_layer = S u / D
T_layer = 2 gamma e / D
R_over_core = R_layer + T_layer^2 R_core / (1-R_layer R_core)
R_patchy = (1-c) R_core + c R_over_core
```

Use expm1 for u, and exact S=0, K=0 and t=0 limits. The denominator contains only
positive terms. An earlier subtraction form failed the nearly conservative/thin
layer energy test; the failure log is retained. No energy clamp fixes that failure.
R_layer+T_layer<=1 and the adding equation is passive. Coating can brighten a
substrate darker than the coat, so darkening is tested where physically implied,
not asserted as a universal law.

F0 and roughness use a passive coated-facet partition with hiding c(1-T^2), and
coating suppresses substrate sparkle. This is an unresolved facet approximation,
not a multilayer specular interference solve. Mineral cleavage/lustre and the
existing wet-film, normal filtering, cavity escape and external weathering-rind
paths are shared. Grain coats are evaluated before cavity trapping and external
rinds; they are not confused with `SurfaceState.coating` on an exposed rock face.
Surface saturation never manufactures a rind.

`without_coatings` is a mass-preserving ablation: it disperses the same coat
mineral stock back into the core mixture. Cement stays partitioned. It does not
remove iron from the body. Reports name this explicitly, so the comparison is a
morphology change rather than a mass-loss experiment.

## Validity, unresolved calibration and visual limits

This is opaque, dry-substrate, three-channel diffuse transport plus the inherited
reduced wet/specular surface response. It omits spectral interference, fluorescence,
polarization, transparent transmission through an entire specimen, explicit 3D
intergrain ray transport, polydisperse phase sizes and calibrated microcrack density.
K and S inferred from RGB are not unique. Inherited l=50 micrometres for many
minerals keeps coarse felsic grains bright. Correcting mixing therefore exposes
large granite/diorite/granodiorite/gabbro/schist overbrightness. The lane deliberately
does not fit new per-rock albedos or replace unknown mineral constants to conceal it.

The inherited v2 material-coordinate grain/band pattern is a bounded zero-mean
optical contrast diagnostic. Its phase selection is not a mineral stock geometry
solver. In sediment patches it does not draw literal cement menisci or individual
micrometre shells. Spheres show mean diffuse plus fixed-light GGX shine; CPU patches
show the shared diffuse/normal structure. Neither is an Unreal screenshot.

Basalt is not exempted from the physical switch. Tholeiite changes slightly, while
the 0.6 mm alkali island-cap realization changes substantially and remains a flagged
bright failure. Marble and obsidian move little because their single/dominant-phase
areal and intimate endpoints nearly agree. Chalk and organic-rich shale remain
intimate and preserve Y to floating precision. All moves are in the per-realization,
three-footprint ledger. An unwelded ignimbrite is explicitly a pale ash scenario;
its high Y is not relabelled as welded rock to pass a target.

Offline reference recollections: Kubelka and Munk (1931), two-flux layer theory;
Hapke (2012), *Theory of Reflectance and Emittance Spectroscopy*, 2nd edition;
Cornell and Schwertmann (2003), *The Iron Oxides*, 2nd edition; Tucker (2001),
*Sedimentary Petrology*, 3rd edition. These support mechanisms, not the chosen
numerical priors. No pages, quotations, web verification, source extraction or
spectral-library measurements are claimed. design review F3/F6 and the prior lane's v3
change ledger provide the local review evidence. Numeric acceptance bands are
unverified sanity ranges, not confidence intervals.

## Reproduction

Use offline Cargo with CARGO_BUILD_JOBS=8 and the lane's persistent target directory.
`tools/build_v4.py` writes only `data/optics_v4.json`. `rockspec_export_v4 DIR`
exports all 67 realization comparisons at 0.1, 2 and 50 m, 134 CPU patches and the
128 migrated diagnostics. `tools/render_v4.py DIR` produces paired realization
rows. The migration retains the v2 surface histories and footprints, takes nominal
modes from explicit v3 reference realizations, and retains the declared altered
basalt mode. Each migrated input and its original v2 input are exported together.
No direct reading of lane-local artifacts is required by library tests.
