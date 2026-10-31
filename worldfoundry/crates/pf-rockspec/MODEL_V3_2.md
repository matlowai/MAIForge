# pf-rockspec v3.2 material structure

Explicit opt-in: construct `v3_2::Evaluator::new(&resolved_body)` once for a
resolved v3 or v3.1b body, then call `evaluator.surface_sample(&surface, footprint)`.
The caller must retain that body and rebuild the evaluator after a body revision.
No environment switch, default-catalogue change, new realization or new material
admission occurs. The sample's evaluator contract is `pf-rockspec-structure/3.2`;
its original MaterialRef still names its admitted v3/v3.1b material. This separates
an optical algorithm revision from geological identity.

This is the native shared material crate used by engine render consumers. The
CPU sheets exercise that actual crate and resolved bodies. They are manufactured
catalogue diagnostics, not an alternate educational app or generated world.
The current Unreal shader still runs v3. A v3.2 shader port and its hardware
parity/visual gates remain consumer work outside this crate's write scope.

The v1, v2, v3, v3.1b and v4 data, algorithms and exporters remain on their original
path. Two optional internal hooks select v3.2 structure preparation and normal
retention. Both are `None` for every pre-existing public evaluator. Neither the
old serialized body nor the old SurfaceSample layout gains fields.

## Inputs and evidence

The source is the resolved body's closed mineral mode, grain distribution,
structure descriptor, fabric, welding grade, independent joint sets, and
fabric-to-material quaternion. Coordinates are metres in the owning body's rest
material frame. Changing a tile, surface ID, camera or fragment renderer does not
change the material seed. Fragment motion must advect the rest frame.

`Record` is the versioned, exported v3.2 structure/fabric record. Its geometric
aspect is inherited from `spec.grain.aspect_ratio`, including the existing
Estimate's reference and plausible support interval. The original v3 `Fabric`
record has spacing and directional strength ratios, **no geometric aspect field**.
Those strength ratios are not converted into an invented shape law. The new
record makes the grain-to-fabric shape proxy explicit.

All added numeric choices are project priors, not extracted measurements.
Each carries value, unit, plausible support, distribution/choice, conditions,
reference and evidence kind. Exact constants in normalization, hashes, integration
and bounds are mathematics, not empirical material parameters. Relevant priors:

| Quantity | Nominal / support | Meaning and evidence |
|---|---|---|
| Foliated cell aspect | inherited, presently 3:1 for schist/phyllite | v3 grain prior, interpreted as in-plane / perpendicular size |
| Cell equivalent diameter | supplied grain d50 for foliated cells | schist 1.5 mm, phyllite 80 micrometres; their old structure wavelengths were 1.5 and 6 mm respectively |
| Fiamme flattening | `1 + 11 W`; support `[1+3W,1+19W]` | grade W in [0,1]; affine process proxy, not a calibrated welding law |
| Band boundary jitter | 0.4; support [0.2,0.45] of nominal spacing | declared bounded stochastic thickness model |
| Inclusion radius | uniform [0.16,0.44] of cell size | mean 0.30; broad marked-size population, not recovered porosity |
| Inclusion occupancy | 0.65; support [0.4,0.8] | independent Bernoulli thinning, optical sites only |
| Outcrop contrast budget | 0.25; support [0.1,0.4] | convex combination with grain-scale deviations |
| Trace half width | 0.08; support [0.04,0.12] of interval width | optical relief support, **not joint aperture** |
| Mica tangent tilt | 0.25; support [0.1,0.5] | uniform per-cell tilt around fabric normal |
| Mica selection fraction | existing model prior 0.12, [0,0.3] | eligible pearly phase cells only; no added specular energy |

RS32-MODEL bibliography, recalled offline and **not consulted or verified in this
lane**: Fossen, *Structural Geology* (2nd ed., 2016), for fabric and joints;
Passchier and Trouw, *Microtectonics* (2nd ed., 2005), for oriented platy grains and
metamorphic fabrics; McPhie, Doyle and Allen, *Volcanic Textures* (1993), for
flattened pumice/fiamme and vesicular textures. These support qualitative model
choices. None supplies the numerical intervals above. Exact chapter/page claims,
specimen calibration and literature extraction are deliberately absent. Existing
catalogue reference IDs retain their original offline qualifications.

## Grains and foliation

Use the inverse fabric quaternion to transform position and face normal into the
fabric frame. Schistosity/crenulation select foliated cells even when an old
structure enum says Laminated, as it does for phyllite. Micaceous also selects
them. The cell axes for equivalent diameter d and aspect A are
`[d A^(1/3), d A^(1/3), d A^(-2/3)]`. The shape is oblate, volume preserving, and
its short axis is the supplied fabric normal. Phyllite uses its actual 80 um d50,
not its inherited 6 mm lamina wavelength.

Full-cell jittered Voronoi sites select mineral optical phases by supplied volume
fractions `area_i = (mass_i/rho_i) / sum(mass/rho)`. The crystal search visits at
most 125 cells, offsets -2 through +2. It is exact: an own-cell site is at most
sqrt(3) away; every omitted cell is at least 2 away. Distance to a candidate cell
box prunes impossible candidates without removing possible winners. Exact ties
use the lower site key. There are no mutable random generators or size-biased
nearest-neighbour weights. Non-foliated crystalline/clastic cells use the old
structure wavelength as their nominal scale. This remains an optical tessellation,
not discrete mineral stock or a liberation simulation.

## Aperiodic bands and fiamme

Band boundaries at nominal scale d are `b_i = i + U_i[-0.4,0.4]`, with a hashed
body phase shift. Adjacent widths have a triangular distribution on [0.2,1.8],
mean 1 and variance `2 * 0.4^2 / 3`. Adjacent widths are correlated because they
share a boundary; this is a displaced-boundary process, not independent renewal
sampling or a fitted lognormal bed model. A floor and two comparisons locate an
interval, including negative coordinates. No scan over all prior bands is needed.

Each interval independently selects a mineral phase with the same volume weights.
Width and phase use different hash salts. Repeated adjacent phases are allowed;
there is no mandatory light/dark sequence. A smooth, hashed 2D value field bends
the micro-band coordinate by 0.3 d over 6 d support. Translating z by nominal d
neither repeats the phase nor fixes the thickness. The warp is a declared optical
fabric irregularity; it is not a deformation-history solver.

Welded ignimbrite bypasses infinite phase bands. Existing glass-bearing modes
supply the lens/matrix optical contrast. A body welding grade W gives oblate axes
with long/short ratio `1+11W` and unit volume determinant, scaled by the supplied
structure wavelength. Compact 3D ellipsoidal kernels give finite lens-shaped
sections with lateral terminations and a size distribution. The nominal welded
record W=1 gives 12:1. Unwelded W=0 does not silently receive welded fiamme.
These are correlated optical domains inside the existing material. No new glass,
pumice mass, void volume or depositional event is created.

## Pores without narrow-jitter lattice rows

The optical inclusion process uses one potential site per full unit cell, a
Bernoulli occupancy mark, a position uniform over the **whole** cell, and an
independent radius uniform on [a,b]=[0.16,0.44]. All 27 neighbouring cells are
visited. Spheres can cross their placement-cell boundaries. Their compact kernels
are summed, not clamped or max-combined:

`P(x) = sum occupancy * max(0, 1 - |x-center|^2 / radius^2)^2`.

This is a thinned, stratified marked point process, not exact Poisson or hard-core
packing. It removes the strong lattice spectral line of v2's narrowly jittered
centres. Residual short-range stratification is a stated model limit. Because all
radii are below 1/2, at most eight cells contribute at any point, so `0 <= P <= 8`.
This conservative bound keeps overlapping pores passive without saturation or
hidden mean-changing clipping. It also limits pore albedo contrast; normal relief
carries much of their visible shape.

For uniform centers and marks the ensemble/large-volume mean is analytic:

`mu = 0.65 * (32 pi / 105) * (b^4-a^4) / (4(b-a))`.

The gradient is the analytic derivative of the same kernels. A vesicular surface
uses `mean + amplitude * filter * (mu-P)`, with amplitude bounded against [0,1]
and proportional to the supplied vesicularity. Fiamme use the same compact
geometry with oblate axes and actual glass-minus-areal optical contrast, centred
with `P-mu`. Affine volume preservation leaves the mean unchanged. The porosity
and cavity transport already in v3 are unchanged; this is **not** a second void
fraction or a new geometry/collision representation.

## Outcrop term and filtering

A separate set of fields uses actual `fabric.resolved_spacing_m` and each actual
joint-set spacing and plane normal. Fabric is a weak plane, never an opened joint.
Columnar fabric has no single planar fabric trace; its declared joint sets supply
the optical traces. Fabric strength is `1-tensile_ratio`; joint strength uses the
supplied persistence fraction `persistence/(persistence+spacing)`. Plane/face
crossing multiplies the strength by `sqrt(1-dot(n,axis)^2)`, so a parallel plane
cannot paint a false intersection trace. Persistence controls the closure's
support, not finite fracture geometry.

Each irregular interval has a compact quartic trace centred at t=1/2 with half
width w=0.08. Its average over **every complete interval**, regardless of that
interval's thickness, is `mu_trace = 16 w / 15`. The signed albedo term is
`mu_trace - trace(t)`. Normal relief comes from its analytic gradient and the
existing half-fracture-slope prior. The projected trace is integrated analytically
over the scalar pixel footprint, with an explicit point limit for widths below
1e-5 interval units. Six fixed interval visits suffice: within the active compact
filter, a footprint spans less than two nominal intervals and can intersect at
most four jittered intervals. This avoids narrow-trace pixel aliasing without
requiring a stochastic ray bundle.

The old compact filter `max(0,1-(footprint/(2 scale))^2)^2` applies independently
at the micro minor-axis scale and each outcrop spacing. Micro filtering remains
a reduced approximation, not an exact 2D pixel integral. Fine phyllite grains and
millimetre fabrics correctly disappear when unresolved; metre joint fields can
remain on a 10 m patch whose pixels are 10/128 m wide. A 10 m **patch** and a
10 m **pixel footprint** are different quantities. There is no invented coarse
mineral tint to keep unresolved crystals visible.

The substrate reflectance is a convex combination of bounded, centred micro and
outcrop deviations around the existing v3 bulk mean (75/25 when outcrop terms
exist). Each outcrop term receives at most 1/N of its contrast budget. Thus
passivity follows from the amplitude bounds rather than post-hoc clipping.
Far footprints set every deviation and resolved normal to zero exactly. The v3
cavity, unresolved roughness, F0, coating and wetness closures remain the bulk
limit. Thin-film normal retention is applied separately at the micro wavelength
and each outcrop spacing; a submillimetre film must not erase metre relief.

The mean claim concerns linear substrate reflectance, ensemble/large-volume
sampling and the exact far-footprint closure. A finite face patch, nonlinear wet
reflectance transform or lit image need not equal the bulk mean. This is not a
white-furnace proof of a full multiscale BRDF, a field calibration or a mechanics
model. The existing reduced roughness closure is retained, not replaced by a new
measured slope spectrum. Exported finite-patch means make sampling error visible.

## Mica through a standard one-normal renderer

On a pearly mineral cell, a separate material hash accepts a facet with probability
`model.sparkle_fraction`. A deterministic, bounded tilt around the fabric normal
supplies its orientation. Orient the facet into the face hemisphere and project
it into the face tangent plane with slope bounded by 1.5. This cap is a numerical
model choice for cross-sections, not a crystallographic angle measurement. The
same physical micro filter removes the facet when unresolved. Coatings and water
retain or suppress its resolved normal through their own state.

The selected facet replaces that cell's optical normal perturbation. Standard
GGX uses the unchanged F0, roughness and diffuse energy. No emission or extra
highlight energy is added. `structure_v2.flake_specular_fraction` is zero in this
path because the selected normal already represents the facet. Existing aggregate
`Sparkle` metadata remains descriptive; do not use it to add another light lobe.
`flake_normal_material` records the unprojected candidate orientation for diagnosis.
The qualification string records whether a cell was selected, for diagnostic
masks only. A renderer should use the actual normal, not parse that string.
Appearance still depends on lighting and roughness. This is a single-normal
projection, not a polarized/multilobe mica BRDF or a camera-keyed twinkle effect.

## GPU translation contract

`portable_structure()` exports the version, exact seed as low/high u32 words,
derived record, per-phase areal weights and raw single-phase RGB response,
pearly/glass masks, joint sets and quaternion. The existing material transport
and surface-state tail inputs are still required. No position-dependent baked
textures are needed. This is a portable input record, not a GPU parity receipt.

Use the exact v2 integer permutation (multipliers 0xbf58476d1ce4e5b9 and
0x94d049bb133111eb), wrapping u64 arithmetic or the existing uint2 implementation,
and length-delimited FNV over body_id then realization_id. Preserve signed-cell
bits and floor for negative coordinates. Uniform is `(hash >> 11) / 2^53` in the
f64 reference. An f32 shader may round uniforms; boundary-pixel phase flips must
be measured, not ignored. No tile, chart, camera, traversal order, surface ID,
view direction or time enters the seed.

Per-sample bounds are fixed: up to 125 crystal cells, 27 inclusion cells, at most
64 phases, at most 8 joint sets plus one fabric term, and 6 integration intervals
per trace. The normal-retention hook repeats only trace work when a coating or
film is present. The CPU can skip impossible crystal cell boxes. Tables should
be prepared once per material revision; per-phase transport coefficients are
already precomputed by Evaluator. No rejection loop or unbounded renewal process
exists. GPU performance is not measured here.

The f64 reference refuses active scaled coordinates outside `abs(p/scale)<2^40`.
This is not an f32 precision promise. A GPU port must preserve integer cell plus
fractional coordinate, using a stable rest-frame origin or split coordinates.
Do not cast planet-scale metres to f32 or change the material seed on origin shifts.
At near-boundary points, report numerical mismatch with the same common footprint
and coordinate precision used by the reference.

## Validation and diagnostics

Run offline from the repository root, with a persistent target directory:

```sh
CARGO_TARGET_DIR=external-data/target CARGO_BUILD_JOBS=8 cargo test --offline --release --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml -- --nocapture
CARGO_TARGET_DIR=external-data/target CARGO_BUILD_JOBS=8 cargo run --offline --release --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml --bin rockspec_export_v3_2 -- external-data/replay-v3_2
python3 worldfoundry/crates/pf-rockspec/tools/render_v3_2.py external-data/replay-v3_2 unpublished-design-record review-lanes/mat-structure/replay-figures
```

The exporter writes paired v3/v3.2 patches for all 67 v3 realizations at 1 mm,
0.1 m, 2 m and 10 m widths, each 128 pixels square. It exports beauty, base
reflectance and normals, 268 point probes and the portable records. Extra mica
faces use 20 mm schist/granite and 1.6 mm phyllite patches, at 256 pixels square.
Lighting is fixed one-normal GGX plus diffuse, no per-material exposure fit.
Pillow only labels existing CPU images; no network or image model is used.

Tests cover exact nearest-site search against a larger neighbourhood, jittered
width statistics and repeat phases, the pore kernel's analytic mean/gradient and
lattice Fourier component, trace integration against independent quadrature,
all 67 realizations and 8 footprints including wet/coated samples, exact far bulk
identity, frame rotation, reload/fragment identity, changed-body contrast,
foliation anisotropy, finite welded lens terminations, outcrop visibility,
statistical substrate means, invalid inputs and single-normal mica selection.
Old tests and exporters remain separate compatibility gates. Current lane results,
retained failed candidates and source hashes are in the lane's RESULT and INDEX.
