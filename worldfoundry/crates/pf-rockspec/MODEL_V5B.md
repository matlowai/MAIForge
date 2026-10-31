# Optics v5b: explicit coherent and particulate forms

This opt-in evaluator separates intact cut/broken/outcrop surfaces from an explicitly sized particulate surface. It uses the existing engine `v3::ResolvedBody` and the v5 USGS mineral calibration. It never changes admitted mineral stock, identity, density or mechanics. Runtime remains procedural. The existing v5 data, evaluator outputs and version remain frozen.

```rust
let optics = pf_rockspec::v5b::Optics::embedded()?;
let form = pf_rockspec::v5b::SurfaceForm::Coherent {
    finish: pf_rockspec::v5b::SolidFinish::Broken,
    additional_rms_slope: 0.2, // caller's declared dimensionless slope, not a fitted default
};
let optical_view = optics.prepare(&admitted_body, form, Some(persisted_optical_seed))?;
let sample = optical_view.surface_sample(&actual_face_state, footprint_m)?;
```

The result identifies `pf-rockspec-optics/5b` and `usgs-form+connected-chord-km/5b`. The `MaterialRef` still identifies the admitted v3 composition. The engine's face/material consumer must supply form from actual state. No automatic inference from a rock key, camera or desired colour occurs. Cut/broken/outcrop names do not themselves change colour: coherent cases with identical slope and body inputs are identical. No game/world consumer is silently migrated. Structure 3.2c and its GPU-port inputs remain untouched.

`SurfaceForm::Particulate` requires a positive finite particle diameter and an explicit additional slope. It describes an optical counterfactual or an already-admitted particulate material, not a comminution event. Applying it to a solid body does not create fines, change volume or debit stock. A consumer creating actual powder must use the engine's stock and material-state mechanisms separately.

## Why v5 needed a form dimension

V5 infers effective absorption K and scattering S from aggregate mineral reflectance factors. Its scattering is `1/d + 1/l`, with inherited internal scattering lengths. The latter term often dominates at coarse grain sizes, so a coarse solid retains much of the reflectance of a particulate reference. V5b preserves the inferred absorption and changes the density of scattering interfaces according to explicit optical form.

The colour method remains exactly v5: USGS cvASD samples, 380–780 nm, analytic CIE1931 colour matching, illuminant E, Bradford adaptation to D65, linear sRGB. The limitations of that approximation remain. These coefficients are not a uniquely recovered microscopic single-scattering albedo or a calibrated bidirectional reflectance function.

## Coherent-domain closure

For each mineral in the calibrated optical view, recover the existing coefficients:

```
K_i = ((1-R0_i)^2 / (2 R0_i)) * (1/d0_i + 1/l_i)
S_before_i = 1/d_old + 1/l_i
```

`d_old` is the admitted optical grain divided by its existing boundary efficiency. For a coherent crystalline phase, use chord `d_i=d_old` and the retained fraction of powder-style internal scattering

```
h_i = L / (L + d_i)
S_after_i = 1/d_i + h_i/l_i
```

The closure tends to the original floor for fine/disrupted domains and to interface-density scaling proportional to `1/d` for large connected domains. L is an effective scale of the retained interface population, not a measured crack depth or film thickness. This is a reduced-order preparation/coherence closure. It is not a first-principles solution of all grain contacts or an assertion that scattering inside every real crystal follows this law.

For coherent glass, there are no crystalline grain boundaries. The idealized internal term is removed, `h_i=0`, while inherited absorption is retained. Approximate the connected phase chord by the mean run length of occupied cells:

```
d_i = min(d_old/(1-v_i), max(d_old, d_cap))
```

At `v_i=1`, use `max(d_old,d_cap)` directly. The fraction `v_i` comes from admitted dry mass fractions and mineral densities. `d_cap=0.01 m` is a declared optical correlation prior with support `[0.001,0.1] m`; it is neither measured body thickness nor generated geometry. Crack/inclusion scattering and refractive-index matching at real glass/crystal contacts remain unresolved. The cap is never selected from a rock-key brightness target.

For particulate material, cap the physical optical domain by the supplied particle diameter and restore the particulate internal-scattering term:

```
d_opt = min(d_old, d_particle)
h_i = 1
S_after_i = 1/d_opt + 1/l_i
```

Physical grain detail is capped separately by `min(admitted_grain,d_particle)`. The optical pattern wavelength is capped in the local view for particulate rendering. This is not a rewrite of the source grain distribution or body stock.

Per-phase reflectance is `KM(K_i/S_after_i)`. The exact old K is retained across form queries for the same body and seed. The v4 mixing engine still evaluates separate volume-weighted K/S and areal endpoints with v5's frozen transition. Since the intimate/areal weight can change, a mixed material need not obey a universal monotone brightness law even though each phase becomes brighter with more scattering. The lane reports a mixed-basalt exception rather than clipping its result to obtain a form-ordering pass.

## Reusing the engine evaluator without changing legacy paths

The new optical view represents its positive per-phase K/S using an equivalent common bookkeeping diameter:

```
d* = max(2/min_i(S_after_i), 2*d_physical)
l*_i = 1/(S_after_i - 1/d*)
R*_i = KM(K_i/S_after_i)
```

The common diameter is not a new geological grain size. With reference diameter d* and scattering length l*, the unchanged v4 coefficient equations reproduce K_i and S_after_i. Its optical boundary-efficiency field is adjusted only in the local optical view. Receipts expose physical phase chords, original/new scattering and unchanged absorption; they are the meaningful physical records, not the bookkeeping diameter.

The only edit to v5 is making its existing borrowed material-view method crate-visible. Its method body and all numeric operations are unchanged. The other existing source edit registers the v5b module. All old data files are unchanged. The lane regenerated every v5 export and the v4 numeric exports to verify this boundary.

Existing coating/cement stock partitions and actual supplied face coatings still run through the common evaluator. Coating geometry and unknown mineralogy are not invented from a dark reference swatch. Existing partitions on a particulate optical counterfactual are inherited approximations, not a solved redistribution of coats during crushing.

## Roughness and escape

Additional dimensionless RMS slope is explicit, supported in `[0,1]`. It adds unresolved variance to the existing roughness. It also provides a passive diffuse escape approximation:

```
g = 1/(1+slope^2)
R_surface = R*g / (1-(1-g)*R)
roughness_out = min(1, roughness_in^4 + slope^2)^(1/4)
```

This sums repeated cavity returns; it cannot add diffuse energy. Normal-incidence Fresnel coefficients are retained, and diffuse energy is recomputed as `R_surface*(1-F0)`. Layer and structure receipts are updated consistently. It is a microcavity approximation, not a full angle-resolved scattering model. Diagnostic slopes 0, 0.2 and 0.4 are declared cut/broken/outcrop scenarios, not empirical finish constants. The fit uses zero *additional* slope; the engine's existing surface/structure response still operates.

## Calibration and qualification

Only one usable local USGS record has the combination of known fresh solid preparation and a sufficiently relevant catalogue proxy for fitting this closure: fresh basalt interior BR93-46B, record 11848. Its exact modes and roughness are not known, so this remains conditional. White chert has unknown form in the available documentation; the limestone target is weathered. Neither is promoted to an uncoated form-calibration target.

A fixed 41-point search over `L=10^(-5+i/10) m`, i=0..40, selects **0.001584893192461114 m**. The three-channel solid-target root-mean-square error is **0.0044699566**. Three documented USGS physical powder mixtures remain checks of the particulate branch, with RGB errors **0.0319071**, **0.0592326** and **0.0317173**. They do not identify coherent L. The configuration retains the full tested L support as a support interval, explicitly not a confidence interval. No glass absorption or glass-chord cap was fitted in this arc.

Parameters were frozen before v5b independent evaluation. RELAB colours, spread and associated preparation metadata remain lane-local validation; no RELAB value is stored in the model or used to choose a coefficient. The model still misses some powder and glass-rich references. A reduced overshoot does not constitute a completed glass calibration. Refer to the lane's failed gates before considering default adoption.

USGS Data Series 1035 distinguishes particulate mineral measurements from solid rock measurements and their different integrating-sphere geometries. Primary preparation evidence and caveats are recorded in the lane source notes. The fitted retention equation and finite connectivity cap above are explicitly this project's model assumptions, not formulas attributed to USGS.

## Form-indexed distributions

The same optical seed chooses the same v5 mineral observations before applying either form. Thus form comparisons are paired virtual specimens with unchanged composition and absorption. Sampling occurs once during preparation, not once per pixel. It preserves the existing joint RGB observation draws, not independent channel noise.

The lane exports 67 realizations × four form contexts × 64 seeds, with separate means, spreads and covariances. Measured reference populations are retained separately by form **and condition**. Where a reference library describes coated versus uncoated surfaces, those labels remain separate; their measured variation is not converted into extra colour noise or a wider fitted slab factor. These validation-only populations are not runtime lookup tables. The predicted mineral-population spread is not a claim to have measured whole-rock occurrence frequencies or all surface histories.

## Tests and reproduction

The native runner is `rockspec_forms_v5b OUTPUT_DIR` with optional `--no-images`, `--fit-usgs TARGETS_JSON`, or `--requests REQUESTS_JSON`. Requests contain only material realization and explicit form context. They do not contain measured target colours. The fit route accepts the separately audited USGS target packet. The request route never fits parameters.

The default exporter writes 268 procedural patches, 67 realization/form records at three footprints with 64 seeds per form, and 128 diagnostics with original/migrated inputs. All diagnostic face histories and explicit altered modes are retained. Six new tests cover phase K conservation, passive form/roughness response, body and old-output preservation, invalid inputs, 128 migrated histories and the labelled USGS training regression. All 106 crate tests pass.

Review, exact commands, source hashes, every referenced-rock figure and unresolved gates: `unpublished-design-record review-lanes/mat-calibrate/v5b/INDEX.md`. Heavy outputs: `external-data/`. These are diagnostics of the native material subsystem, not a standalone teaching application, stock-transforming physics feature or proof of native-world/GPU adoption.
