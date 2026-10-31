# pf-hydro-catchment

Native routing library for T-Q24 (D-0155). **Runtime integration is pending.**
The retained real C7 check passes joint r8 and coupled r32 head/conservation gates.
The old fixed-parent nested candidate remains a measured failed diagnostic.
The CLI writes evidence and returns an error if the coupled numerical gate fails;
numerical success does not establish visual or streamed-player acceptance.

This crate has its own `[workspace]` and does not change the parent workspace.
It reads pinned C7 products, derives a catchment from any requested planet cell,
and runs one priority flood and metric D8 over its whole r8 mesh. It never
mutates bed elevations, erodes terrain, writes canon, or restricts destruction.
D-0139/D-0141: no geographic allow-list or manually chosen river geometry.

## Domain and boundary definitions

The requested cell's first downstream declared water body (or an explicit closed
floor) identifies a catchment. Full reverse reachability adds **all** upstream
coarse cells, including intervening lakes, then a one-cell D8 halo. Water-body
IDs, stages and outlet metadata come from the pinned endpoint hydrology. Lakes
absorb water at this domain tier; declared outlet edges are retained for future
coupling, not silently extended to the ocean. A receiving lake is not necessarily
endorheic. Other intercepted lakes keep their own identities and budgets.

Fine dry-land beds use the attributed Q20 continuous spherical relief law on
endpoint elevation. The globally seeded relief is a conditional unresolved
initial condition, not evolved channel history. Native wet cells retain their
archived coarse bed and stage. This explicit coarse water-body representation
preserves very shallow C7 lakes that stochastic dry-land relief would erase;
it retains coarse shorelines. No post-routing clipping or bed changes occur.

The outer halo meets explicit external samples whose fixed heads are the maximum
of sampled bed and interpolated published coarse fill. These are **boundary
exports**, not ocean/lake terminals. They require a neighboring-domain solve
before downstream continuity can be claimed. Halo nodes have zero *local tagged
source*: they route this core catchment's water, not complete neighboring runoff.
Source supply includes runoff over every core sample, including direct supply
to native water bodies, in m³/s from spherical area and pinned uniform runoff.

Priority flood computes a nondecreasing-to-upstream routing head. Strictly lower
heads choose the steepest metric neighbor; flats follow flood parents, with no
invented epsilon slope. Uphill **bed** links remain visible in the audit. Filled
heads are an auxiliary drainage surface, not a simulated water-volume/time
solution. Do not draw dry-land channels over high fills as though bed gradients
were downhill or water volumes had been simulated.

## API and products

- `Input::read`, `Hydro::read`: retained C7 data and provenance.
- `Domain::build`: requested cell → complete upstream domain, r8 mesh and halo.
- `route(&Mesh)`: head, receiver, tagged discharge, terminal ownership and audits.
- `nested::solve`: r32 requested cell and an adaptive interface to actual r8
  samples. Head values are copied from the parent, not interpolated onto a
  potentially higher fine bed. Per-port flux and the actual replacement graph
  are checked; parent flood rank is not a physical constraint.
- `coupled::solve`: replace the requested cell with its r32 geometry, retain r8
  heads as lower bounds outside it, then jointly re-solve and reaccumulate all
  domain sources once. Actual transfers correct parent receivers/fluxes instead
  of treating inherited inlet discharges as immutable. On measured C7, all r8
  heads remain exactly unchanged and the old two-node loop is removed.

Unsupported fixed-parent ingress remains a failed diagnostic with its source
volume and exact loop in `metrics.json`. Never splice that diagnostic into r8.
Use the full coupled product and recorded corrections together. The inherited
r8 head floor is an explicit reduced-order assumption: this refinement does not
lower pre-existing outside depression heads when finer terrain exposes an outlet.

`provenance.json` binds world/source hash, catchment, requested cell, resolution,
shape, radius and law version. `inputs.json` contains input paths and SHA256s.
`routing.jsonl` contains every physical node (no display threshold); external
receivers resolve in `boundaries.jsonl`. Fine IDs are south-first raster IDs at
r8. `joint.csv` and `nested.csv` are unsmoothed plotting edges at Q ≥ 10 m³/s.
`coupled.jsonl` contains the full adaptive product, with stable sample addresses
`[resolution, global_fine_id]` including every external boundary sample. Its cache
identity also includes the requested refinement cell. `divergences.csv` records
changed r8 cross-cell receivers; no clipping to coarse receivers is performed.
`metrics.json` records coupled head/receiver/terminal corrections and owns gates.
`nested_rejected.csv` preserves the old fixed-parent diagnostic. Runtime publication is
false pending server integration and review, even if numerical gates pass.

## Run

Use an assigned target/tmp directory and no more than eight build jobs:

```sh
cargo test --offline -j 8 --manifest-path worldfoundry/crates/pf-hydro-catchment/Cargo.toml
cargo build --release --offline -j 8 --manifest-path worldfoundry/crates/pf-hydro-catchment/Cargo.toml
pf-hydro-catchment --c7 C7_ROOT --out NEW_OUTPUT --cell GLOBAL_CELL --refine 8 --nested 32
```

For the retained comparison, `--q20 Q20_RUN` reads its data-selected crop instead
of `--cell`. This is an evaluation selector, never a service allow-list. Run
`scripts/plot.py OUTPUT` with NumPy and Matplotlib; it reads existing Q20 curves
and produces the three-panel crop and requested-cell comparison. No new world
is generated. Real acceptance measurements and commands live in the canonical
river-catchment lane's REPORT.md and COMMANDS.md.

## Attribution and limitations

`input.rs`, grid/relief in `terrain.rs`, and plotting conventions are adapted from
the retained, physically rejected Q20 `pf-hydro-mr` prototype, original source
base `198033663053e43ba3a347dd52d2f2b9943d2b63`. Exact reused file hashes are in the
lane's BASE.json. Q20's per-cell portals, bed incision and material redistribution
are not reused. New code implements routing only. Native/Q20 plotted curves are
398–400 My surface-pass products, while Q24 uses endpoint terrain/hydrology.
Their provenance and geometry differ; this is not an identical-terrain ablation.

Uniform C7 runoff only; no climate expansion, sediment transfers, temporal
meandering, erosion, MFD, lake-volume dynamics or canonical admission. Longitude
wraps; raster rows do not jump across the poles. The domain-size guard is a
resource limit, not a geographic restriction. Streaming adapters must reconcile
this declared bed law with actual served terrain and honor water/boundary status.
