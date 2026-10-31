# Opt-in structure revision 3.2c

2026-10-25. Module `v3_2c`, version `pf-rockspec-structure/3.2c`.
Use `v3_2c::Evaluator::new(&resolved_body)?.surface_sample(&surface, footprint_m)`.
`portable_structure()` exports the immutable coefficients and laws for a GPU
consumer. No default selector changes. `lib.rs` only adds the module and its
implementation of the existing structure-to-optical-tail interface. All older
structure, catalogue, optics and exporter source files are retained unchanged.

This is the shared material subsystem used by MAIForge's rock renderer.
Inputs are the existing resolved geological body, its constituent modes, fabric,
joint sets and material/rest coordinates. An engine consumer must select this
version for the same admitted body. The diagnostic planes are manufactured test
fixtures, not a new world, geological admission or completed gameplay feature.
UE port, measured GPU timing and native in-game acceptance are separate gates.

## Finite traces and cost

The retained physical/conditional inputs are the mean spacing s, finite-disc
diameter L, aperture a, set normal, and two orientation slope supports. L and a
retain exactly the v3.2b supported symmetric uniform distributions. Orientation
retains independent uniform tangent slopes in `[-tan(theta),tan(theta)]`.
The normal is normalized `(slope_u,slope_v,1)` in the supplied set basis.
Fabric coherence remains a declared optical prior: 16s, supplied support [8s,32s],
sampled by the existing symmetric-support rule. It does not terminate actual
penetrative fabric or create physical discontinuities.

Let `A = pi E[L²] E[cos(theta)] / 4`, where `E[cos]` uses the same 16x16
midpoint quadrature as v3.2b. Desired site density is `rho=1/(s A)`. The mean
number of finite-disc intersections per unit scanline normal to the set is 1/s.
Let `rmax=Lmax/2` and `sinmax=t/sqrt(1+t²)` for
`t=sqrt(tan(theta_u)²+tan(theta_v)²)`. Cell dimensions are now

```
dx = dy = rmax + 1.08 s sinmax
dz = rmax sinmax + 1.08 s
lambda_cell = rho dx dy dz
p = lambda_cell / 16
N ~ Binomial(16,p)
```

Reject bodies exceeding this bounded capacity. One uniform draw and a stored
17-entry cumulative distribution function (CDF) produce N with at most five
binary-search comparisons. Then only N sites are generated, each uniformly
jittered over the full cell with independent length, aperture and orientation
marks. Unlike a clipped Poisson draw, this distribution preserves its declared
mean exactly. Its variance is `lambda_cell*(1-lambda_cell/16)`; it is **not an
exact Poisson process**. v3.2b also used a bounded binomial cell count, implemented
as 16 separate occupancy trials. Changing cell dimensions changes pair/count
correlations. Marginal mark distributions and projected-area intensity remain
the same; exact realizations do not. The scanline and distribution tests measure
that distinction instead of claiming identical point processes.

The maximum active pixel footprint is 2s. A patch's bounds in set coordinates
are no larger than the above cell dimensions. Search at most 27 neighboring
cells. At actual footprint f, reject a cell before hashing if its box lies
outside `[rmax+w sinmax, rmax+w sinmax, rmax sinmax+w]`, with `w=.08s+f/2`.
This rejection is exact. It is tested against a 125-cell search and finite kernel
support/gradient checks. The new cells are tied to maximum persistence instead
of generally oversized mean-persistence tangential cells.

Trace hashing is a declared uint32 model revision. `hash32` uses xor shifts
16/15/16 and wrapping multiplies `0x7feb352d`, `0x846ca68b`. Cell coordinates are
mixed with `0x9e3779b9`, `0x85ebca6b`, `0xc2b2ae35`; the existing once/body
64-bit FNV seed is folded by XOR of its low/high words. Uniforms use the high
24 bits divided by 2^24. Count takes the cell hash directly. Each site needs one
hash for its key, three for position, one for length, and (after sphere rejection)
two for orientation and one for aperture. Thus:

| Query work per active set | v3.2b | v3.2c |
|---|---|---|
| Cell candidates | 27 | at most 27 |
| Cell key hash | 4 uint64 hashes/cell | 1 uint32 hash per surviving box |
| Count/occupancy | 32 uint64 hashes/cell | one CDF lookup, no additional hash |
| Sites | 16 occupancy trials; generate accepted sites | directly generate N <=16 |
| Extra hashes/site | 4 before sphere, 3 after | 5 before sphere, 3 after |
| Worst hash calls/set | 3996 uint64 | 3483 uint32 |
| Typical hashes/set, measured 15 sets | see lane cost JSON | about 204 to 576 uint32 at f/s=.02,.5,1.9 |

These are algorithm operations, **not GPU milliseconds**. The old shader's
uint64 emulation cost makes replacing its hashes especially relevant, but only
consumer timing can establish the GPU gain. Crystal and pore/lens functions
retain their 64-bit hashes; costs above exclude them. Coated/wet samples may
query normal-bearing traces a second time through the retained optical tail.
Compositional bands do not repeat a normal query because they add no relief.

The compact LOD (level of detail) multiplier is
`max(0,1-(f/(2s))²)²`. At f>=2s return an all-zero probe before cell/hash work.
This cutoff already existed in v3.2b; v3.2c preserves it and explicitly tests
zero operation counters at and beyond the cutoff. It is not claimed as a newly
invented saving. The zero-integral compensated kernel and its normal-footprint
box filter are unchanged. LOD multiplies a zero-mean field, so removing a set
does not change its bulk optical mean. This closure is not an exact 2D pixel
area convolution, and no unbiased pixel-normal distribution is claimed.

GPU implementation must preserve wrapping uint32 arithmetic and comparison
order, bind the precomputed CDF, precompute tangent slopes/bounds per body,
and keep bounded loops. CPU trace cell coordinates are restricted to absolute
values below 2^30. The GPU must split exact material cell origins from local
float fractions; shifting world positions without preserving material cell
identities changes the field. No GPU parity is claimed for this crate revision.

## Fiamme without columns

The retained lens/matrix partition uses only existing glass and residual body
constituents. It conserves each bulk constituent volume exactly. Welding still
sets the nominal flattening `1+11W` and core fraction `.8W`. Occupancy .65 and
uniform radius [.16,.44] are unchanged project priors, not measured geology.

Host slab boundaries are `i+U[-.4,.4]`, with independent boundary hashes. Width
is triangular on [.2,1.8], mean one. Each slab has an independent uniform in-plane
rotation and full-cell shift. Within it, each occupied XY cell owns one confined
lens. Z is normalized to slab width, with the chain-rule gradient divided by that
width. Lenses cannot overlap within cells or across slabs. There is no shared
vertical column lattice. These spatial choices are serialized in the new record.

For a slab of width h, lens volume is h times the old lens volume; volume of the
slab is also multiplied by h. Width and occupancy/radius are independently keyed.
Therefore the previous mean `mu=.65 J(core) E[r³]` is preserved, as are the
partition and its constituent balance. A particular finite cut need not equal
this mean. Local minor dimensions now vary by [.2,1.8] about the nominal value.
The footprint fade uses the nominal minor dimension. It is an approximate
antialiasing closure; exact lens pixel integration is not claimed.

Tests retain the old lens-layout function as a falsifier. Three seeds integrate
column density through 128 slabs, with 2048 horizontal samples; additional
horizontal and vertical line spectra check that the change does not simply
move the periodicity. Mean/nonzero-variance checks prevent erasing the lenses
to pass. The independent volume/gradient test checks 131072 points and 768
finite-difference derivatives. The old failure evidence stays in lane outputs.

## Bedding, compositional bands and relief

The supplied quartz sandstone and unwelded ignimbrite descriptors already say
`Bedding`, with resolved spacing .12 m, but their tensile ratio is one. Optical
bedding does not require mechanically weak planes. 3.2b's `1-tensile_ratio`
therefore erased valid supplied fabric. New optical strengths are .35 for
`Bedding` and .4 for `FlowBanding`, conditional on a declared resolved fabric.
They are project priors with plausible interval [.1,.7], dimensionless, stored
as Estimate records. No catalogue fabric or mechanical ratio is modified.

Two narrowly read primary sources support the **occurrence**, not these numeric
strengths. Hamblin's radiographic study reports internal stratification even in
apparently homogeneous sandstone, and variable visible expression controlled by
composition, weathering and diagenesis. This supports separating optical
visibility from mechanical anisotropy. [Hamblin 1965, KGS Bulletin 175(1)](https://www.kgs.ku.edu/Publications/Bulletins/175_1/).
Millward describes weak bedding, massive units, and welded eutaxitic textures
within ignimbrite sequences. Bedding is therefore conditional on the body
fabric, not imposed universally on ignimbrite. [Millward 2003, BGS RR/02/02,
Froswick and Lincomb Tarns formations](https://webapps.bgs.ac.uk/Memoirs/docs/B06513.html).
Neither source calibrates our amplitudes, spacings, random process or colours.

Bedding also receives bounded compositional strata at the descriptor spacing,
so nonzero optical traces are not its sole, very faint, cue. Irregular warped
boundaries select existing constituent colours with the same body volume
probabilities, using independent width and phase hashes. Subtract the areal
colour expectation, then limit amplitude with a **body/footprint-dependent**
headroom bound. No selected-pixel-dependent normalization, tint, new constituent
or mineral stock is introduced. This is a conditional optical sorting model,
not an explicit simulation of deposition. Welded `FlowBanding` keeps its fiamme
representation; an extra strong compositional stripe field is not added to it.

The old gneiss band-depth ramp is removed. Gneiss/bedding/flow fabric traces add
no groove normals because differential weathering relief is not provided by
those descriptors. Joint groove normals and existing cleavage/foliation optical
weakness relief remain. Mica still uses one bounded facet normal, not an added
specular lobe. No new weathering geometry is asserted. Tests isolate gneiss from
joints: 1024 normals remain exactly the input normal while colour varies.

## Independent micro contrast and mean preservation

Crystal selection, phases and the original microscopic amplitude remain.
Compute deterministic minimum/maximum colours allowed by that micro field at
this body and footprint. Additional bedding contrast uses at most 75% of the
available headroom. Trace albedo uses 25% of remaining headroom divided by
`27*16*set_count`, a conservative overlap bound. **Do not multiply the micro
field by .75.** Near grain-scale views retain their full crystal mosaic contrast.
The pessimistic trace-albedo bound is intentionally faint; joints mainly read
through their groove normals. Improving that bound requires a separate proof,
not overlap clipping or a fitted mean correction.

Each field is centered independently. Phase/width independence preserves band
means; analytic inclusion volume preserves fiamme means; compensation preserves
finite-trace means. Headroom depends on body and footprint, not chosen phase or
trace value, so it does not reweight the stochastic fields. All surface values
remain bounded without clipping. Once every structure filter is zero, the
retained cavity/variance optical closure and normal are exactly the legacy
result. Finite lit patches and nonlinear coating/water responses need not have
the statistical substrate mean.

## Validation and reproducibility

The lane report provides source hashes, full logs, operation counters and paired
images at spans 1 mm, .1 m, 2 m and 10 m (128 pixels each). Tests require no
`_local` data. To optionally save statistical JSON, set
`PF_ROCKSPEC_V32C_EVIDENCE` to an evidence directory; unset means only that export
is skipped, not the tests. Run `cargo test --offline --release` with this crate's
manifest. `rockspec_export_v3_2c OUTPUT_DIR` renders the same 10 body fixtures,
light, coordinates and optical tail for v3.2b and v3.2c. The accompanying
`tools/render_v3_2c.py` only assembles and labels those CPU pixels with Pillow.

Limits: no source spectrum calibration, no GPU timing/parity, no new mechanical
fractures, no displacement or stock, no fabricated geological identity, no claim
of exact Poisson gaps, no global guarantee of absence of all spectral peaks,
and no default admission. Consumer review must include GPU parity and timing,
wet/coated behavior, scale-safe normals and a native game view before adoption.
