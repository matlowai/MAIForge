# pf-rockspec v3.2b: finite traces and readable fiamme

Opt in with `v3_2b::Evaluator::new(&resolved_body)`, then `surface_sample`.
The new evaluator ID is `pf-rockspec-structure/3.2b`. The admitted material's
original v3/v3.1b MaterialRef is preserved. There is no default switch or new
catalogue. All v1/v2/v3/v3.1b/v4 and v3.2 evaluators remain separate frozen paths.
The only shared implementation change is a private structure trait in lib.rs
that dispatches preparation and layer-normal retention to the chosen revision.
The original v3_2.rs, its exporter, tests, data and label script are untouched.

This remains the native engine's shared material crate. The CPU diagnostics use
its real catalogue and resolved-body interfaces. They are manufactured material
fixtures, not a standalone game/course or admitted generated world. A consumer
must port and qualify this revision before claiming new Unreal appearance.

## Why the v3.2 outcrop fields failed visual review

A displaced one-dimensional boundary process still supplied one infinite trace
per nominal interval. Jitter changed exact widths but preserved a recognizable
rhythm and uninterrupted crossings of joint sets. Persistence only attenuated
contrast; it never ended a trace. Every set used one direction. The integrator's
reported vertical rules and sandstone plaid followed directly from that model.

V3.2b replaces only those outcrop fields with finite, marked plane patches.
The successful v3.2 grain tessellation, micro-banding and pore kernels are retained.
Fiamme are the other explicitly changed structure, described below.

## Declared joint and fabric distributions

`Record.outcrop_distributions` exports each `TraceDistribution`. It copies the
body's spacing, persistence and aperture Estimate records, including units,
references and uncertainty. Spatial variability is declared separately from the
original epistemic prior intervals. An uncertainty interval is not silently
claimed to be a measured spacing distribution.

| Quantity | Spatial model | Source and qualification |
|---|---|---|
| Mean line-intersection spacing s | finite-disc point intensity, targeting 1/s intersections per metre along the nominal set normal | supplied `joint.spacing_m`; conditional model, not measured gaps |
| Persistence P | uniform `[P-r,P+r]`, `r=min(0.75P,P-low,high-P)` | supplied persistence mean/support; finite disc diameter, all values remain inside the prior |
| Aperture a | same supported symmetric-uniform construction | supplied mean/support; fabric stays exactly zero |
| Set dispersion | independent uniform tangent slopes, bounded by tan of half the supplied angular support span | azimuth support projected by sin(polar); polar support gives the other tangent direction; at poles the polar span supplies both directions |
| Fabric optical coherence | mean 16s, support [8s,32s], sampled as the supported uniform law | explicitly new optical prior; no v3 fabric persistence exists |
| Fabric dispersion | tangent half-width 5 degrees, plausible support [0,10] degrees | explicitly new optical prior around the supplied fabric frame |
| Optical relief half width | `s * (0.04 + 0.04*a/(a+0.0001s))` | bounded reduced closure linking actual aperture to trace width; this is not the physical gap width |

The angular half-widths also have typed Estimate metadata. Tangent-u is the
azimuth direction away from poles; tangent-v is the polar direction with a sign
that is immaterial for a symmetric distribution. These are bounded tangent-slope
distributions, not a uniform distribution in Euler angles or a calibrated Fisher
orientation distribution.

Every optical site is a finite circular patch. Its radius is half the sampled
persistence. Its normal is the set normal plus the sampled tangent slopes,
normalized. A quartic lateral taper vanishes smoothly at the patch boundary.
A closed penetrative fabric does not become a collection of open fractures:
its finite **optical coherence patches** terminate, not the geological fabric.
Columnar fabric still uses its supplied joint sets, not a fictitious planar
columnar-fabric trace.

A set uses an anisotropic placement grid in its supplied material frame. Each
cell has 16 independent Bernoulli candidate sites, with centers uniform over the
entire cell. Empty and multiply occupied cells are allowed. There is no one-plane-
per-interval requirement, fixed phase ordering or periodic boundary sequence.
This is a thinned stratified point process. It approximates broad, near-exponential
intersection gaps; it is not advertised as exact Poisson sampling.

For projected mean patch area A, choose site intensity `lambda=1/(s A)`.
For uniform diameter L and dispersed normal, `A=pi E[L^2] E[cos(theta)]/4`.
The cosine expectation is evaluated once with a fixed 16x16 midpoint quadrature.
The per-trial presence probability is `cell_volume*lambda/16`. Invalid capacity
is an error. All supplied catalogue sets fit. The approximation in the cosine
quadrature is separate from the exact mean-reflectance argument below.

Cell dimensions enclose the largest patch projection plus the largest active
pixel footprint. Tangent cells are at least P wide and expand when necessary.
The normal cell width is `s + max_radius*sin(max_tilt) + 0.08s`.
A fixed 27-cell query then includes every possible contributing site. A larger
125-cell oracle checks this bound. The first prototype refused salt/gypsum because
its tangent cell was 15 mm too small; the grid construction was corrected rather
than dropping those rocks or truncating their patches.

## Zero-mean finite trace and normal relief

For normal distance z, half-width w and normalized distance t=z/w, use

`h(t) = (1-t^2)^2 (1-7t^2)` for |t|<1, zero otherwise.

It has a dark center and compensating shoulders, range [-32/49,1], zero integral,
and zero value/derivative at its endpoints. Its antiderivative is
`t - 3t^3 + 3t^5 - t^7`. A box-filtered h and its derivative are computed directly
from endpoint evaluations; very small footprints use the analytic point limit.
The footprint is projected across the actual dispersed patch plane.

Multiply h by `(1-r^2/R^2)^2` for in-plane radius r<R, zero beyond R. The
normal integral is zero for every lateral location. Consequently every finite
patch, at every position, orientation, aperture and persistence, has zero volume
integral. Summing overlapping patches preserves that result. No max-combination,
saturation, overlap clipping or fitted constant is used to obtain the mean.

The original compact scale filter acts at each set's spacing. Plane/face crossing
attenuates a patch parallel to the cut. Contrast also uses the declared persistence
support fraction, or the fabric tensile anisotropy, as in the prior model.

The color budget uses a conservative bound of 27*16 possible overlaps per set,
then the existing 25 percent outcrop share. This makes the albedo component small.
Most visible traces come from analytic **normal relief**, not painted dark
lines. Normal relief uses a separate single negative-height groove
`-depth*(1-t^2)^2` with the same finite lateral envelope and analytic footprint
filter. Its slope is the positive gradient of the unsigned compact kernel.
The color-compensation shoulders are not copied into height: doing so made the
first candidate look like parallel scratches around each joint. The groove has
one minimum and no subsidiary ridges; its compact gradient integrates to zero.
The optical height profile is not published as actual collision geometry. Each joint set contributes its own existing half-fracture-slope amplitude.
Independent sets add relief; they are not weakened by dividing their normals by
the number of sets. Normal vectors remain in the visible hemisphere after tangent
projection and normalization. Coatings and water retain detail independently at
each spacing, as in v3.2's scale-aware normal hook.

Far footprints set all new fields exactly to zero. The original v3 bulk transport,
cavity, F0, roughness, coating and wetness results remain the exact limit. The
claim is about linear substrate means and the far closure. Individual finite
patches, nonlinear wet response and lit images need not have the bulk mean.
This is not a white-furnace proof of a complete multiscale BRDF. Optical patch
positions are not an authoritative fracture network, collision mesh or stock map.

## Fiamme diagnosis and correction

The supplied welded ignimbrite contains about 57.7164 percent glass by solid
volume. Its existing pure-glass single-phase linear RGB is approximately
[0.31803,0.34498,0.35917]. Its crystals are brighter. Thus retained mineral optics
already permit contrast. The weak v3.2 image was not evidence of identical glass
and matrix optics.

V3.2 nevertheless treated fiamme as additive soft pore kernels. The average
kernel weight was only 0.02046, with a conservative maximum of eight overlapping
contributions. That bound reduced allowed contrast, and the rounded kernels had
little homogeneous interior. Their thin dimensions also limited resolved area.

V3.2b instead declares material lens domains. Each lens stays within its host
cell, so supports do not overlap and the exact field range is [0,1]. Radius and
occupancy retain v3.2's declared priors. The in-plane placement orientation and
global phase are hashed from the material address, never the camera or tile.
The fiamme's physical oblate axes retain the welding-driven ratio `1+11W`.
For W=1 that remains 12:1. No tint or new mineral is introduced.

A uniform glass core has radius fraction `c=0.8W`. A cubic smoothstep in squared
radius joins this core to the outer lens boundary. The volume integral is analytic:

`J(c)=4pi [c^3/3 + integral_c^1 r^2 * (3u^2-2u^3) dr]`,
where `u=(1-r^2)/(1-c^2)`.

The mean lens fraction is `mu=occupancy*J(c)*E[radius^3]`. At W=1, mu is
0.0667149653. The supplied existing glass species, in their original relative
proportions, make up the lens mode. The residual matrix's constituent volumes are
`matrix_i=(bulk_i-mu*lens_i)/(1-mu)`. Insufficient existing glass is an error. For the supplied
body the reconstructed original phase volumes close exactly in f64 in the test.
The matrix and lens dry-mass modes follow from those same densities. These are
conditional optical domain fractions, not newly admitted discrete mineral stock.

A profile/partition-only prototype raised 0.1 m luminance standard deviation from
0.00607658 to 0.00780276, a factor of 1.284, below the preset 1.5 gate. Its evidence
and images are retained. The remaining limitation was local connected-glass
transport: both domains still used the same body-averaged boundary efficiency.

The added local healing closure is explicitly parameterized by composition and
welding, not display color:

`boundary_survival = 1 - W*(1 - bulk_glass_volume_fraction)`

`lens_optical_diameter = existing_bulk_optical_diameter / boundary_survival`.

Interpretation: as welding completes, a pure-glass fiamme core removes the
crystal/shard interfaces represented in the bulk-averaged boundary count. This is
a declared reduced transport prior, not a measured welding-rate law. The interval
[bulk_glass_volume_fraction,1], units, conditions and evidence kind are exported.
The matrix keeps its existing optical diameter. Mineral n/k, scattering priors,
RGB tables, the body mode and the body bulk mean are unchanged.

For the supplied W=1 body, survival is 0.5771641 and the optical diameter changes
from 20 mm to about 34.65 mm in the lens. The same existing mixture function gives
lens RGB [0.279743,0.305670,0.319447], versus residual matrix
[0.407292,0.414462,0.409934]. There is no hand-assigned color or gain multiplier.
The unchanged structure contrast prior is applied to the lens-minus-matrix
response and centered with the analytic mu around the retained bulk mean.
This centering preserves the legacy bulk closure; it is not a new heterogeneous
radiative-transfer solve.

A welding 0/0.5/1 control checks that the exported coefficient changes through
this formula and that the ordinary mixture function produces the exact lens
response. An identical-glass lens/matrix control must produce exactly zero optical
contrast. This would fail for a tint applied irrespective of composition.

## Evidence status and assumptions

All new quantitative closures are explicit project priors. The old input
Estimate records remain unchanged. RS32B-MODEL refers to this document and the
retained MODEL_V3_2.md bibliography. Those books were recalled offline, not
consulted in this lane. No web search, literature extraction, measured joint
statistics, specimen calibration or new empirical confidence intervals are claimed.

In particular, finite circular patch shapes, supported-uniform marks, tangent
orientation laws, fabric coherence, optical relief width, homogeneous lens cores
and the local welding/connection law are conditional abstractions. They remain
subject to geological and native-render review. The v3.2 bulk-transport limits
and lack of a complete physically calibrated BRDF also remain.

## Reproducibility and GPU bounds

The f64 reference uses the unchanged exact length-delimited FNV material seed and
64-bit permutation. No tile, surface ID, camera, time or traversal order is hashed.
Material-space frame and negative-coordinate floor rules are unchanged.
`portable_structure()` exports every new distribution and lens coefficient.

Per active trace set, work is bounded by 27 cells times 16 candidate sites. There
are at most eight joint sets and one fabric term. Most candidates are rejected by
occupancy and a bounding sphere. There is no rejection sampling, renewal scan,
recursive search or unbounded loop. The inherited micro crystal and pore bounds
remain 125 and 27 cells. Layer-normal retention repeats trace evaluation only when
coating or film is present. GPU cost and f32 parity are not measured here.

Preserve integer cell and fractional coordinates when porting. The active f64
scaled-coordinate domain remains below 2^40; this does not justify raw planet-scale
f32 positions. A native port must measure phase-boundary and normal discrepancies
at the same footprint and coordinate precision before enabling this version.

The periodicity falsifier uses 2048 samples over 128 nominal spacings for eight
rocks and three body addresses. The largest non-DC Fourier-bin power divided by
all positive-frequency power must be below 0.10. Field standard deviation must
exceed 0.005, preventing a pass by erasing the field. A periodic sinusoid must fail.
Finite-patch output normal tests independently require nonzero actual render
normals at the 10/128 m footprint. These are stated engineering falsifiers, not a
claim to have proved the absence of every possible visible artifact.

Run all tests and the comparison exporter offline:

```sh
CARGO_TARGET_DIR=external-data/target CARGO_BUILD_JOBS=8 cargo test --offline --release --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml -- --nocapture
CARGO_TARGET_DIR=external-data/target CARGO_BUILD_JOBS=8 cargo run --offline --release --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml --bin rockspec_export_v3_2b -- external-data/review-v3_2b
python3 worldfoundry/crates/pf-rockspec/tools/render_v3_2b.py external-data/review-v3_2b unpublished-design-record review-lanes/mat-structure/v3_2b/review-figures
```

The exporter compares v3.2 and v3.2b for ten affected/control bodies at 1 mm,
0.1 m, 2 m and 10 m patch widths, using the same 128x128 positions and footprint.
It produces beauty, base reflectance and normal images, portable inputs, point
probes and supplementary mica controls. It does not change the old exporter.
The lane's INDEX/FINDINGS contain source hashes, failed-candidate evidence,
resource measurements and exact legacy/v3.2 export comparisons.
