# pf-rockspec

Opt-in material structure v3.2 is documented in [MODEL_V3_2.md](MODEL_V3_2.md).
Use `v3_2::Evaluator::new(&resolved_body)` for foliated grains, irregular bands,
fiamme, jittered pores, outcrop traces and single-normal mica facets. Existing
evaluators and data stay on their original paths. The v3.2 Unreal port remains
a separate consumer gate.

**V2 is opt-in:** see [MODEL_V2.md](MODEL_V2.md) for connected-domain optics,
family/scenario-specific coatings, unresolved joint/cavity transport and
material-addressed mesostructure. V1 remains the default and is preserved.
The following original documentation describes **v1**; its old limitations
are historical v1 evidence, not a description of all v2 behavior.

Standalone, low-memory **candidate shared engine material catalogue**. Own Cargo
workspace and lockfile. No parent-workspace membership or dependencies on other
MAIForge crates. No feature is enabled in the world or Unreal by this crate.
The integrator must bind admitted material realizations to the actual world.

There are 43 exact shared rock keys and 32 mineral/constituent records. The latter
include the requested groups plus epidote, talc, apatite, cordierite and an explicit
low-absorption felsic glass endmember. `organic_matter` is a constituent, not a
mineral species. Clay, amphibole, pyroxene, serpentine and glass records are group
priors, not one chemically exact species. No new shared **rock** keys were added.

Every table quantity has a value, interval, unit, interval kind, distribution,
conditions, reference, evidence kind, support scale and covariance group. All are
**project priors**, not measured specimen values or source-extracted endpoints.
Bibliographic references motivate the properties and relationships. The exact
interval endpoints were chosen for this implementation and require qualification.
No network or source downloads were allowed. Do not cite these ranges as measured
literature confidence intervals. Correlation labels identify joint constraints;
there is no fitted covariance matrix or independent prior sampler.

## Engine contract

`GeologicalMaterialSpec` supplies mineral mass modes and ranges, clast/crystal
size percentiles, optical domain scale, crystallinity, intrinsic porosity,
interparticle porosity, vesicularity, density, mechanical priors, fabric, joint
sets, and weathering pathways. Thermal laws and weathering rates are explicitly
unknown. Mechanical fields use `Property::KnownPrior` or `Property::Unknown`.
Unconsolidated deposits have no fictitious intact UCS, tensile strength, E,
Poisson ratio or fracture toughness. Their dry friction angle is a broad prior.
Mechanics are dry intact-specimen priors, not a failure law; directional failure,
confinement, moisture dependence and calibrated strengths remain consumer work.

`MaterialRef` contains spec ID/version, realization ID and body ID.
`Catalogue::realize_prior` is an **explicit diagnostic/conditional realization**,
not automatic material admission. The caller retains the resulting `BulkState`,
reference and realized optical grain scale once per material domain. Production
must replace conditional modes with source-supported or admitted priors. The
sampler never rerolls composition, infers ore, samples a random seed or reads a
renderer colour preset. No fracture/weathering consumer is wired yet.

`BulkState` stores the current dry solid mineral mass fractions (including any
buried alteration), intrinsic pores, deposit packing pores, saturation,
temperature, damage, a unit fabric quaternion, and alteration provenance.
`SurfaceState` separately stores surface and creation-event IDs, optional elapsed
exposure time, forcing record, rest-material coordinate, unit face normal,
fracture/cleavage/abrasion mechanism, wet coverage, water-film thickness,
abrasion, roughness history, and an optional explicitly supplied coating layer.
An unknown elapsed time remains `None`. Age alone has no optical effect. A fresh
cut retains the supplied altered bulk. Existing exterior/joint faces retain
whichever coating and history their caller supplies.

`surface_sample(&material, &surface_state, footprint_m)` returns:

- Linear diffuse base reflectance RGB, energy-partitioned diffuse coefficient,
  Fresnel F0 RGB, and perceptual GGX roughness (`alpha = roughness²`). These are
  pre-lighting quantities, not display sRGB.
- Separate substrate and coating optical samples and coating weight. A coating
  never mutates bulk mode or density.
- Substrate and coating normal-detail bands at their **separate physical
  wavelengths**, oriented in the material frame. Unresolved slope variance goes
  into roughness; the model does not invent a blended physical grain diameter.
- Optional substrate and coating sparkle budgets. Each is a redistribution of
  specular energy, **not extra radiance to add on top of F0**.
- Material/surface IDs, model version, intrinsic dry bulk and packed-deposit dry
  density, and the conditional-prior qualification.

No time, environment, external file, global mutable state, renderer or texture
is consulted by `surface_sample`. Output is deterministic for pinned input on
the same numerical platform. Cross-architecture bit-identical `exp`/`sqrt` is
not claimed. Future GPU parity uses the supplied outputs and the FINAL_PLAN
numerical tolerances, not a claim that a backend already passed.

## Physical conventions and reduced models

Mineral modes are **dry solid mass fractions**, not QAPF volume proportions or
surface coverage. Mass fractions close to one; unknown/malformed modes fail
rather than normalizing away lost stock. Specific solid volume is `sum(w_i/rho_i)`.
For a statistically random section, area/volume fraction is proportional to
`w_i/rho_i`; this weights the optical mixture. Density is
`rho_dry_intrinsic_bulk = (1-phi_intrinsic)/sum(w_i/rho_i)` and packed-deposit
density is that value times `(1-phi_interparticle)`. Vesicles are a **subset** of
intrinsic pores, not an additional density subtraction. Resolved cavities are
excluded upstream. Water remains a separate stock. Optical evaluation changes
neither stock nor mineral identity.

Cleavage directions/angles are mineral crystallography. Rock joint sets describe
spacing, persistence, aperture and normal polar/azimuth angles relative to the
fabric frame. They are broad conditional input descriptors for a fracture lane,
not generated planes or site observations. Cleavage quality encodes weaker and
polymorph-dependent sets. Listed angles are representative distinct angles, not
an exhaustive pairwise matrix. Foliation is a rock fabric, not a mineral species.

Mineral diffuse RGB is an optically thick **clean granular aggregate** prior at
1 mm optical domain scale, not the reflectance of one polished transparent
crystal. An isolated quartz cleavage/fracture face reflects approximately its
Fresnel value; the much larger diffuse coefficient here comes from repeated
internal/particle scattering. Effective n is an orientation/composition average;
plagioclase also records An mole-fraction support. Group-specific n, density and
An are correlated priors and must be resolved together. Opaque mineral records
also carry broad extinction-index k priors. Their n/k values are not recovered
visible spectra or calibrated ore reflectometry.

The optical equations are:

1. Invert diffuse half-space reflectance to absorption/scattering ratio,
   `q = K/S = (1-R)²/(2R)` (Kubelka-Munk two-flux approximation).
2. A scattering closure `S(d) proportional to 1/d + 1/l_internal` scales q by
   `(1/d_reference + 1/l_internal)/(1/d + 1/l_internal)`. Grain size increases
   absorption path until internal scattering dominates. These lengths are
   explicitly uncalibrated priors. Clast size in conglomerate, gravel and till
   is separate from the smaller internal optical domain.
3. Mix q by constituent volume fraction. Recover
   `R = 1/(1+q+sqrt(q²+2q))`. This is stable and monotonic. It assumes a
   homogeneous intimate scattering mixture, not resolved granite grains or
   gneiss bands. It may over-darken mixed rocks with small opaque fractions.
4. Compute constituent Fresnel `F0=((n-1)²+k²)/((n+1)²+k²)` and area-average.
   For dielectric k is zero. Roughness comes from declared lustre priors and
   separate fracture, cleavage, glass, abrasion and fabric slope terms.
5. A footprint attenuates a grain-scale normal band by
   `d/sqrt(d²+footprint²)`. Resolved and unresolved squared slopes sum to the
   original variance. Unresolved variance contributes to GGX alpha². This is
   a declared analytic low-pass closure, not a measured spatial power spectrum.
   Mean reflectance does not change merely because the camera footprint changes.
6. For an explicit coating, weight is
   `coverage * (1-exp(-2*extinction*thickness))`. Blend separately evaluated
   substrate/layer coefficients. This is a rough absorbing-layer/coverage
   approximation, not a full multilayer adding-doubling optical solve. A zero
   thickness/coverage layer is transparent. No lithology ID is blended.
7. For wet areas a Lambertian return has trapped fraction `p=1-1/n_water²`.
   Summing repeated returns gives `Rwet=R(1-p)/(1-pR)`, following the internal
   reflection mechanism described by Lekner and Dorf. This is a declared
   simplified model, not a full implementation of their paper. Explicit surface
   wet coverage mixes dry/wet areas; pore saturation alone does not tint a face.
   A water film progressively replaces the top Fresnel lobe and suppresses
   micro-slope variance over a calibrated-in-future smoothing length.
8. Report diffuse energy `(1-F0)*R`. F0 plus that coefficient cannot exceed one.
   Full hemispherical BRDF energy conservation is not established by this
   coefficient check. Thin-film interference, detailed wet-interface refraction,
   temperature optics and polarization are not modeled.

The initial single dark volcanic-glass surrogate made pumice/ash implausibly dark.
The final table declares a separate `felsic_glass_low_absorption` endmember in
felsic volcanic mixtures, with explicit prior uncertainty. Dark `volcanic_glass`
is an inclusion/iron-rich optical surrogate appropriate to the chosen obsidian
and mafic examples. A world generator must supply actual glass composition; it
must not choose endmembers just to make a desired colour. There is no universal
albedo for either glass group. Porosity affects density but **vesicle-wall optical
transport is not solved**. The pale pumice example is therefore unqualified as a
quantitative model of pumice, even though its mean appearance improves.

Optical validity is an optically thick, non-emissive, clean/rough/coated mineral
surface near ambient conditions and visible RGB. Transparent thin halite/calcite
crystals, hot lava, subsurface transmission, birefringence, biological crusts,
fracture geometry, cavity shadows, layered gneiss, grain-resolved colour masks,
hydration/perlite microcracking, and time-varying reaction laws are unsupported.
Temperature, saturation and damage remain state inputs for consumers; no optical
law for those fields is silently invented. Inputs are validated numerically, but
caller evidence policy must also enforce this physical validity domain.

## Data and reproducibility

- `data/catalog_v1.json`: versioned catalogue and full reference qualifications.
- `data/shared_keys_v1.json`: the exact 43-key interface.
- `data/diagnostic_inputs_v1.json`: 128 frozen manufactured inputs, with full
  material refs, bulk states, histories, normals, layer data and footprints.
- `data/swatch_inputs_v1.json`: 43 rocks x 3 histories x 3 filter footprints.
- `tools/build_catalog.py` and `tools/build_diagnostics.py`: retained deterministic
  provenance generators. Editing them requires versioning/review before adoption.
- `rockspec_export`: CPU-generated output JSON, validation and sanity tables.
- `tools/render_swatches.py`: consumes those outputs only. It does not reimplement
  the material evaluator. Existing local Pillow is used; Rust needs no image crate.

```sh
export CARGO_BUILD_JOBS=8
export CARGO_TARGET_DIR=external-data/target
cargo test --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml --offline
cargo run --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml --offline \
  --bin rockspec_export -- unpublished-design-record review-lanes/mat-rockspec/evidence
python3 worldfoundry/crates/pf-rockspec/tools/render_swatches.py \
  unpublished-design-record review-lanes/mat-rockspec/evidence
```

The contact sheet uses one neutral light and identical smooth ellipsoids. Labels
show rock ID, dominant mass-mode minerals, linear base luminance and roughness.
The 0.1 m, 2 m and 50 m columns are **averaged evaluator footprints**, not detailed
photographs of hand samples/outcrops/cliffs. At these scales most mineral detail
is already unresolved. The same mean colour across scale is intentional. No
random grain speckle, texture or macro fracture geometry is invented for the
picture. Real cleavage morphology and cliff-scale variation require downstream
geometry and spatial material realization. The sheet alone cannot qualify those.

The three approximate sanity bands are supplied in the operator brief (basalt
0.05-0.10, granite 0.25-0.35, white marble >=0.6), not independently verified
literature measurements. The marble upper bound 0.9 is a declared loose project
sanity bound. Matching those checks is necessary evidence, not empirical optical
validation. The lane report records misses, visual defects and exact resources.
