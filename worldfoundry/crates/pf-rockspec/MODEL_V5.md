# Opt-in USGS optical calibration v5

`v5::Optics::embedded()?.prepare(&resolved_v3_body, seed)?.surface_sample(&face, footprint)` evaluates the shared engine material substrate with a calibrated optical view. `seed=None` uses mean per-mineral reflectance; `Some(u64)` draws one complete RGB observation per mineral once for the body. Persist and inherit that seed with the body's optical realization if a consumer opts in. Reprepare after a bulk-state change. Sampling never changes identity, mineral stocks, density, mechanics or geological placement. No spectral file, texture, library lookup, or random draw occurs per pixel. This is not an automatic world/Unreal migration.

The output contract is `pf-rockspec-optics/5`, numerical model `usgs-rgb-km+grain-transport/5`. The supplied body remains `pf-rockspec/3`. Existing source evaluators, data and exports remain unchanged. Only the new module is registered in lib.rs. The 3.2c structure contract and GPU-port inputs are untouched. v5 reuses v4's shared surface and coating evaluator; the diagnostic exporter uses its existing v2 structure, not a newly selected 3.2c structure or a visual redesign.

## Evidence and colour method

All shipped fitted numbers derive from USGS splib07, public domain, Kokaly et al. 2017, DOI 10.5066/F7RR1WDJ. Each observation stores the original record, specimen/fraction ID, title, relative numeric file and SHA256. RELAB contributes no fitted number and no runtime distribution. It is independent lane-local validation only.

The fixed colour method integrates USGS splib07b cvASD 1 nm samples over 380–780 nm under equal-energy illuminant E. Colour matching uses the Wyman, Sloan and Shirley (2013, JCGT 2(2)) analytic CIE1931 2-degree approximation. Bradford adaptation converts E to D65, followed by XYZ to linear sRGB. There is no display transfer curve in model coefficients. Diagnostic pictures alone apply sRGB encoding. This reuses the supplied reader's method, with full precision rather than its four-decimal output. It is an explicit approximation, not certified tabular CIE colourimetry. Linear-sRGB luminance can differ slightly from the original illuminant-E Y quoted by the integrator.

Missing/sentinel, negative or above-one visible samples are rejected, never interpolated across missing visible channels or clipped for fitting. Relative-reflectance and Nicolet-only observations are rejected. Mineral fits accept purity a/b/c, retaining weak-contamination spread as uncertainty; d/unknown purity is excluded. Only the mineral chapter is used for mineral fits, excluding named mixtures and films. Duplicate instruments for a specimen/fraction are not independent draws: prefer usable Beckman, then a usable ASD measurement. The admission count is 41/56 minerals, 374 specimen/fraction observations, with 65 numeric-size titles. Fifteen minerals keep their resolved-v3 priors.

Title size ranges use arithmetic midpoints and upper limits use half the limit. These are labelled nuisance assumptions, not measured medians. Explicit nominal sizes are used directly. HS fraction codes remain undecoded. A spectrum without a numeric size retains the old reference-diameter prior. Composition, purity, preparation and grain differences all contribute to empirical spread; its min/max and covariance are not confidence intervals.

## Per-grain inference

For an observed channel reflectance factor R, infer the effective two-flux ratio

```
q = K/S = (1-R)^2 / (2R)
q_at_d0 = q * (1/d_observed + 1/l_prior) / (1/d0 + 1/l_prior)
R_at_d0 = KM(q_at_d0)
```

For unknown size, d_observed=d0. For each mineral, average these transferred RGB vectors, or choose one joint vector with a deterministic FNV-1a/SplitMix-finalizer hash of the caller's seed and mineral ID. All three channels use the same observation. This preserves observed hue covariance. The optical catalogue is cloned only during preparation, not during each pixel evaluation. Per-mineral receipts retain record IDs, transferred means, min/max support, covariance and known-size counts. Unknown covariance is explicitly labelled unknown; its zero storage sentinel is not evidence of zero uncertainty.

This is calibrated *aggregate per-grain reflectance*. It is not a unique microscopic single-scattering inversion: reflectance-factor geometry and unknown scattering lengths cannot identify that quantity. The inherited grain-transfer law remains a reduced-order prior. Its extrapolation beyond measured sizes is not validated. No monotone grain-size law is asserted for all natural opaque/specular specimens.

Already-admitted feldspar conditions restrict eligible records. Plagioclase uses the named albite/oligoclase/andesine/labradorite/bytownite/anorthite interval containing the body's An fraction. K-feldspar uses a declared optical bin `(R-G)/Y > 0.10` for pink and its complement for white. This threshold is a modelling condition, not an inferred iron content or a measured universal pink classification. If no observation matches, retain the resolved-v3 prior and report that gap.

## Mixing calibration and retained coatings

V4's separate volume-weighted K and S, areal endpoint, penetration-length transition, mass-conserving coating partitions and finite passive slabs remain the optical mechanisms. The fitted v5 coefficients use transition limits `[0.25,32.25]` and a fine-oxide effective scattering cross-section factor of `8` for hematite and goethite.

The factor scales *both* K and S at the prepared body's optical grain. Thus it changes an oxide's contribution to an intimate mixture without changing its pure-mineral reflectance. Given original S=1/d+1/l, represent S'=f S with l'=1/(f S - 1/d), rebase the reference diameter to d, and preserve R(d). It is an effective particle cross-section correction conditional on this grain representation, not an iron-oxide absorption tint or a measured mean free path. Tests exercise pure-grain invariance at three diameters.

The fit evaluates 324 candidate parameter triples: transition start 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64; width 1, 2, 4, 8, 16, 32; oxide factor 1, 2, 4, 8, 16, 32. Three physically prepared USGS mixtures receive weight 1; fresh basalt and white chert receive weight 0.25 because catalogue modes are not specimen-matched. Documented mixture endmembers are selected explicitly for the fit. A numeric-size nuisance prior is still needed for the ground mixture and for ground goethite WS222. Calculated AMX areal spectra are excluded as independent physical observations. The five-target weighted RGB RMSE is 0.04809, compared with 0.15828 for the same calibrated/matched inputs using factor 1 and v4 transition limits. White chert retains 0.10014 RGB RMSE. Using the default mineral populations instead of the documented endmembers gives 0.08270 weighted RGB RMSE, including 0.13463 for the hematite/quartz mixture. Halving and doubling the assumed powder sizes gives 0.05075 and 0.06034 respectively at the fixed selected parameters.

Transport is weakly identified: 41 grid candidates lie within 5% of the minimum, with transition starts spanning 0.25–64 and ends 5–96. All those candidates use oxide factor 8. The minimum lies on a grid boundary; these limits are not established universal material constants. No held-out performance is implied by the training residual. The lane includes leave-one-target-out grid sensitivity and an independent RELAB comparison.

USGS limestone CU02-11A is a weathered surface according to its original metadata. Its earlier fresh-proxy fit was rejected transparently; old evidence remains in the lane. It is only a condition-labelled comparison now.

All six v4 shell/cement morphology profiles remain priors. Available coated-rock spectra do not provide matched substrate, coverage and thickness measurements. A coating-thickness fit is therefore unidentifiable here and is **not claimed**. The oxide/quartz laboratory mixtures constrain intimate scattering, not shell thickness. Bright-rock and glass-rich rock errors remain. Do not promote v5 as universally calibrated rock albedo or make it the default based on the test count.

## Distributions and consumers

`Optics::calibration()` exposes mineral observations, their empirical statistics, and a complete 49-key rock reference ledger. Each populated rock record carries USGS observations and covariance/spread; unsupported keys explicitly say that no rock-specific fit is available. Full source JSON also retains condition and title metadata. No measured rock mean is applied as a tint. Uniform specimen/fraction weights describe this curated library, not geological occurrence frequencies. Reference-family coverage does not prove the subtype of all 67 realizations.

The lane distribution export supplies every realization's available USGS reference groups plus 64 deterministic propagated mineral draws. The default evaluates mean mineral parameters; nonlinear mixing means its result need not equal the mean of sampled rock outputs. The latter is a model distribution, not measured whole-rock population variance and not a confidence band. Missing-mineral prior uncertainty is not silently inferred from those draws. Minerals with no compatible spectrum retain their prior value, so an all-unreferenced glass body has no measured colour-spread claim.

## Reproduction

Build offline with `CARGO_BUILD_JOBS=8` and a persistent lane `CARGO_TARGET_DIR`. No additional dependencies are needed.

1. `python3 -I tools/build_v5.py COVERAGE_DIR USGS_CVASD_DIR LANE_DIR OUTPUT_JSON [PINNED_FIT_JSON]` extracts USGS data. With no fit argument it emits the unfit factor-1 starting data and five fit targets. With the pinned fit it emits the shipped calibration and its receipt. No RELAB reader exists in this tool.
2. `rockspec_fit_v5 TARGETS_USGS_JSON OUT_JSON` scans the conditional model; optional `--selected` replays the shipped parameters and `--population-check` audits generic mineral populations without refitting; numerical targets and source metadata remain external evidence supplied explicitly.
3. `rockspec_export_v5 OUTPUT_DIR [--no-images]` emits 67 realizations at three footprints, 64 optical seeds per realization, 201 paired CPU patches when enabled, and all 128 diagnostics. The exported original inputs and migrated inputs preserve face history and the explicit altered-basalt mode.
4. `cargo test --offline` runs all 100 crate tests, including five v5 tests. Tests use embedded data, never implicit lane-local files.

The lane report, exact commands, fit evidence, source hashes, coverage gaps, validation comparisons and figures live at `unpublished-design-record review-lanes/mat-calibrate/INDEX.md`. Model source is inside the engine material crate. Diagnostics are evidence of that evaluator, not a standalone teaching application or proof of native-world consumer adoption.
