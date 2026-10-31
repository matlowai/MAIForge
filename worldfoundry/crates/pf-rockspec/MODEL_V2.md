# Version 2 material optics and structure

Opt-in: `Catalogue::embedded_v2()`, `DiagnosticSet::frozen_v2()` and
`DiagnosticSet::swatches_v2()`. `embedded()`, `frozen()` and `swatches()` still
select v1. Package 0.2.0 supports catalogue `pf-rockspec/1` and `/2`. V1 source
JSON is retained byte for byte. Its evaluator path and JSON outputs remain the
key-off reference. New optional fields are omitted from v1 serialization.
No parent crate, engine adapter or world generation is changed.

## What changed and what the parameters mean

1. **Connected-domain optics.** V1 applies the scattering rate of disconnected,
   clean powder at every rock's crystal diameter. That is inappropriate for
   intergrown chert and compact fine mud rocks. V2 declares a boundary-scattering
   efficiency `e` in (0,1], separate from visible grain dimensions. The effective
   transport diameter is `d_optical/e`. The v1 two-flux closure then uses
   `S proportional to e/d_optical + 1/l_internal`. Internal scattering still
   bounds path lengths. No darkening is caused by a rock key in the evaluator.
   `e=1` recovers the disconnected-boundary model. The chosen efficiencies are
   uncalibrated, with explicit prior intervals; they are not specimen measurements.
   A connected domain does not enlarge the mineral grains or change mass.

2. **Conditional alternative modes.** The v2 catalogue explicitly adds minor
   organic, Fe oxide and mafic phases to selected sediment, sedimentary and
   quartzite examples and a magnetite accessory to rhyolite. The alternative
   gneiss is less mafic than v1. Mode sums and derived density are recalculated;
   the evaluator cannot alter them. These are new hypothetical specimens, NOT
   a discovery of those constituents in a world source. V1 shale and mudstone
   already had 5% organic matter. Their correction is principally transport,
   not an invented claim that carbon was missing. Mechanical priors are retained
   as uncalibrated family envelopes, not recomputed from appearance.

3. **Unresolved geometric trapping.** A joint's plane orientation determines its
   crossing factor `1-(n_face dot n_joint)^2`; its persistence and spacing give
   support `p/(p+s)`. A declared symmetric-face RMS relief closure uses half the
   existing fracture/glass slope. The unresolved variance for that joint is
   `crossing * support * relief^2 * f^2/(f^2+s^2)`, with an aperture/spacing term.
   The half-slope is a common model approximation, not measured joint roughness.
   Grain relief uses aspect ratio, packing porosity and d90. Variance enters GGX
   alpha-squared, without changing stock or grain size. The diffuse escape factor
   is `g=(1-0.75*vesicularity)/(1+variance)`. The cavity assumption is that a deep
   vesicle has one quarter opening solid angle; this is not a universal vesicle
   solution. Multiple wall returns give `R_effective=R*g/[1-(1-g)*R]`. That
   includes interreflection and its absorption, and has exact flat/black/white
   limits. Coatings share that geometry. Extra vesicle and pore wall slope
   variance makes vesicular glass matte rather than polished glass. Thin water
   films retain this macro variance and joint variance; they smooth microfacets.

   This is an azimuth-averaged effective transport closure, NOT the complete
   Hapke shadowing law or a fitted directional BRDF. Explicitly rendered joints
   must not also receive the unresolved term for those same scales. A future
   renderer needs a declared resolved/unresolved geometry split. A sample with
   identical geometry descriptors cannot darken merely because its name changes.

4. **Material-frame mesostructure.** The seed hashes length-delimited body and
   realization IDs. Positions and normals use the inverse fabric quaternion.
   Tile ID, chart, surface ID, cut event, time, camera and traversal order are not
   random seeds. Jittered 3D cells represent crystal/clast support. Cell sizes
   vary with the declared d90/d10 spread (bounded radius approximation). Laminae
   and bands are ordered constituent intervals normal to the fabric, with small
   continuous warping. The upper-decile population supplies porphyritic coarse
   structure. Explicit descriptor wavelengths supply band, lamina and vesicle
   spacing. Fine grains are not enlarged to make a distant image prettier.

   Mineral optical contrast is centered around the volume-weighted areal mean
   and bounded so no channel leaves [0,1]. This zero-mean contrast is added to
   the connected-KM bulk response. It is a **reduced heterogeneity closure**, not
   a solution of radiative transfer between resolved neighbouring minerals.
   Expected constituent coverage follows the supplied mass-to-volume conversion;
   a finite field has sampling error. The compact filter
   `max(0,1-(footprint/(2*wavelength))^2)^2` removes it exactly at coarse scale.
   Tests independently sample the near field and bound its mean error, and check
   exact 50 m homogenization. The filter is an approximation, not an exact box
   integral or a proven spatial reconstruction filter.

   Optical grain/clast cells are a projection of admitted modes, not new ore,
   collectible objects or individual geological truth. This arc does not invent
   lithic clast-source identities. Optical normal relief does not displace a
   mesh, alter a silhouette, generate rock joints or supply collision. For
   vesicles random compact spherical pore profiles supply concave wall normals and resolved trapping contrast, centered using their analytic volume integral; finite planar patches have sampling error; full pore geometry is still absent.
   Micaceous patches redistribute a fraction of existing F0 into a narrow flake
   lobe, attenuated by coating/water coverage. No extra additive specular energy.

## Exposure is supplied state

`v2::diagnostic_weathering` requires a typed `ExposureScenario`, v2 material and
an existing face. It creates a **manufactured diagnostic** residue state with
separate composition, thickness, coverage and event/forcing description:

- Wet temperate: mafic Fe oxide residues; partial kaolinite/illite bleaching on
  feldspathic families; impure grey carbonate solution residue; organic-rich
  slaking residue; dark organic oxidation on coal. Soluble evaporites have no
  automatically invented coating. These are example states, not a rate model.
- Arid varnish: explicitly supplied external Mn/Fe/clay dust, with pyrolusite as
  a new constituent. Coal retains its dark organic residue. External deposition
  requires separate forcing and mass accounting before production admission.

The sampler never invokes this helper from elapsed age. Fresh cuts preserve
buried alteration and discard a rind only when the caller supplies a fresh face.
The substrate reflectance is unchanged by a supplied coating; abrasion and
roughness history may independently change its roughness. Zero thickness makes
that coating optically transparent, not a reset of other history.

## Evidence, scientific comparison and limitations

The v2 exporter writes a 2x2 ablation of old/new modes and powder/connected
transport for every rock, with v2 geometry held fixed. This diagnoses mechanisms;
it cannot establish empirical calibration. Every added table number retains the
value/unit/conditions/interval/reference/evidence schema from v1. New endpoints
are project priors. No Bayesian posterior or measured covariance is claimed.

Primary numerical comparison: Lee, Jeong & Kim (2021), *Analysis and Comparison
of Rock Spectroscopic Information Using Drone-Based Hyperspectral Sensor*,
Journal of Engineering Geology 31, 479-492,
https://doi.org/10.9720/kseg.2021.4.479,
https://www.engeojournal.org/articles/xml/b2nd/ , Table 3. The report tool records
all ten measured specimen RGB-band means. Its cited spans are across wavelength
bands, not uncertainty or rock-population ranges. Its basalt and sandstone are
not identified as our precise shared subtypes. No numeric source band is invented
for the other 33 rocks. Basalt and mudstone specimens disagree materially with
integrator expectations; those disagreements are retained in the all-rock table.

`base_reflectance_linear` is a reduced pre-light diffuse coefficient. Rec.709 Y
is not spectral hemispherical albedo, a field photograph, nor the arithmetic
spectral mean in that study. A pass against a rough band is only a sanity result.
A source-span miss is not evidence of failed physics by itself.

## Reproduce

Use `CARGO_TARGET_DIR=external-data/target`,
`CARGO_BUILD_JOBS=8` and cargo `--offline` throughout.

- `python3 tools/build_v2.py` rebuilds only `catalog_v2.json` from immutable v1.
- `cargo run --release --offline --bin rockspec_export_v2 -- OUTPUT` exports
  128 diagnostics, 387 swatches, both scenarios on all 43 rocks, ablations and
  129 material-plane PPM images. `--detail` uses 1024 pixels for each 2 m patch
  (heavy outputs belong on the assigned weights volume). `--no-images` gives the small deterministic sets.
  Exported input sets are the reproducibility reference for the checked-in v2
  diagnostic/swatch inputs. No tests depend on a `_local` path.
- `cargo run --offline --bin rockspec_export -- OUTPUT` regenerates v1 separately.
- `python3 tools/report_v2.py OUTPUT` writes the 43-row JSON and Markdown table.
- `python3 tools/render_v2.py OUTPUT` uses the retained v1 renderer, exported CPU
  outputs and locally installed Pillow. It writes nine comparison pages, the
  full v1-beside-v2 sheet and three closeup pages. No raster-generation service.
- `cargo test --offline` and `cargo clippy --offline --all-targets -- -D warnings`.

The default 0.1/2 m texture images use a physical pixel footprint of span/128;
`--detail` uses span/1024 at 2 m. Their
50 m panel is a uniform, homogenized 50 m cell, not a detailed image with 0.39 m
pixels. The sphere column labels are evaluator footprints. This distinction is
necessary: swatch width and pixel footprint cannot be silently interchanged.
No standalone teaching app, world feature, Unreal backend or native A-D
acceptance is delivered by these crate diagnostics.
