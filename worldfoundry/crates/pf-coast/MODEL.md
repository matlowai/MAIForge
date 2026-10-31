# pf-coast v2: connected coastal potential, opt-in prototype

This is a native T1b consumer in World Foundry. It does not evolve sediment,
modify terrain, create vegetation, or certify an actual delta. Grain size,
coastal rock strength, sediment flux, tidal harmonics and nearshore bathymetry
are absent. The eight-class result is a **potential classification**. The
literature gives useful regime boundaries, not a universal inverse solution
for these eight classes from the available inputs. Full physical acceptance is
therefore open. No threshold was fitted to C7 pictures or class shares.

## Inputs and identity

`coast=1 windows=SET t0=RUN climate=CLIMATE out=NEW tide_range_m=2` explicitly
activates the CLI. Missing `coast`, or `coast=0`, exits before any input read or
output write. No existing producer calls this crate. No pf-region bin keys or
pf-serve spans change. Tide range is required; 2 m is the declared assumption
used in this lane, not a published C7 measurement. Wind is read from the
climate's `wind_u_m_s` and `wind_v_m_s`, never replaced by a constant. The annual
model wind is treated as sustained 10 m wind, an explicit unresolved height,
duration and climatology approximation. The native sea datum is read from T0.

Each immutable context is a set of published windows on the same chart, seed
and recipe. Every consumed field's BLAKE3 is checked. Native fields consumed:
`height_mm`, `fixed`, `shore_mm`, `exposed_litho`, `drainage_area_m2`.
No slope field is published in these inputs; derive the maximum D8 gradient
on the 32 m lattice at each land post touching marine water. This includes
height above the sea datum at the interface. It is not an offshore bed slope.
The input halo is 64 posts; the output halo is 16 posts. Input terrain was
evaluated with the existing 8.192 km margin. No terrain falsifier is changed.

Overlap conflicts in height, sea/lake identity, shore distance, lithology or
**channel** area fail. Area on free posts is not consumed: it is neither a
mouth nor a sediment flux. This distinction is essential: site 5 has 207
free-post area disagreements, but zero channel-area disagreements. No input
bytes are repaired, clipped or averaged. Lake code 4 is excluded, sea code 2
is distinct, only channel codes 1/5 supply mouths. Climate and terrain digests
must match the window upstream list and their grids/seeds must agree.

Context membership is part of every output's identity. Every output window is
evaluated independently from the same immutable context. Exact overlaps hold
within that context. Independently chosen contexts can change fetch, so a
streaming consumer must not mix them. This is not yet an unbounded global
coast service. No invisible mosaic winner resolves conflicting consumed data.

## Laws and primary sources

1. Wind sea: USACE **Shore Protection Manual (1984), eqs. 3-39 and 3-40**, deep
   water limit, with SI wind-stress factor `UA=0.71 U10^1.23`:
   `X=gF/UA^2`, `Hs=0.283 UA^2/g tanh(0.00565 sqrt(X))`,
   `T=7.54 UA/g tanh(0.0379 X^(1/3))`.
   [USACE original manual, hosted copy](https://luk.staff.ugm.ac.id/USACE/USACE-ShoreProtectionManual2.pdf).
   Fetch is traced along 16 rays from a land/marine interface. Upwind rays
   receive squared-cosine weights; wave energies are averaged, periods are
   energy weighted. This directional closure is a project approximation, not
   an SPM prescription. Ray samples are 32 m to 2.048 km, then 256 m to the
   declared 128 km budget cap. Within the context, only `fixed=2` transmits
   marine waves. Outside it, bilinear global T0 elevation below sea is a
   **coarse topographic water proxy**, flagged. It may include inland water.
   The spherical stereographic inverse is extended analytically for these
   global samples without requesting tiles outside the native chart domain.
   Cap hits are separately flagged. No swell, storms, refraction, shoaling,
   depth-limited breaking, duration limit or bottom friction is solved.

2. Relative tide range: `RTR=TR/Hs`; sheltered flat eligibility requires
   `RTR>15`. Masselink and Short (1993), *JCR* 9, 785-800, define RTR with
   **breaking** wave height and distinguish flats above 15. Using wind-sea Hs
   as breaker H is an explicitly approximate regime proxy because marine
   bed depths are absent. We do not implement their fall-velocity beach-state
   model without grain size.
   [Author-supplied paper](https://www.researchgate.net/publication/279908459_The_effect_of_tide_range_on_beach_morphodynamics_and_morphology_a_conceptual_beach_model).

3. Morphological slope boundaries: rocky cliff at >=30 degrees; flat
   eligibility <=3 degrees; gravel-compatible slope >=10 degrees. These use
   NOAA ESI shoreline descriptions. They describe typical morphology, not a
   grain-size measurement or proof of active wave erosion.
   [NOAA ESI v3, sections on slope and ranks 1, 2, 7, 9](https://response.restoration.noaa.gov/sites/default/files/ESI_Guidelines_V3_2002.pdf),
   [ESI v4, gravel rank 6, 10-20 degrees](https://response.restoration.noaa.gov/sites/default/files/ESI_Guidelines.pdf).
   The algorithm's single 10-degree sand/gravel split and application of the
   low-gradient boundary to delta potential are declared simplifications.
   A steep soft bank is included under the requested cliff code; the code
   does not assert impermeable rock where the native lithology says sediment.

4. Width: simple horizontal intersection `W=(TR+R2)/tan(beta)`. On slopes
   below 10 degrees use Stockdon et al. (2006):
   `L0=gT^2/(2*pi)`,
   `R2=1.1 [0.35 beta sqrt(HL0)+0.5 sqrt(HL0(0.563 beta^2+0.004))]`.
   [USGS paper record](https://pubs.usgs.gov/publication/70030520).
   SPM significant period is used as a proxy for the peak period needed by
   Stockdon, with no measured spectral conversion available. This adds
   uncertainty. This empirical runup expression is for natural sand beaches. Its use for
   eligibility widths is an approximation. On steeper banks use H as the
   stated wave-excursion scale, not a claimed runup law. Width is horizontal
   intertidal-plus-runup envelope, **not deposited sediment thickness or
   predicted beach stock**. Above 512 m, width is censored and flagged.
   Zero slope produces a censored finite width. Cliff width is zero.

5. Rock erodibility: native LithologyV1 K multipliers are retained in a local
   explicit table: `[.5,.4,.8,.6,.5,2,.8,.9,3,1.5,3,.7]` for codes 0..11.
   These are the existing model's project priors, **not literature-calibrated
   coastal abrasion coefficients**. Physical motivation for strength-dependent
   erodibility is Sklar and Dietrich's abrasion work, not a calibration of
   these numbers. [Sklar and Dietrich 2004](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003WR002496).

6. Supply: declared channels touching marine posts identify potential mouths.
   Their published area is **effective drainage area Q/r_ref**, not necessarily
   geometric watershed area. A dimensionless proxy is `S=K*A/(K*A+Aref)`, with
   `Aref=1 km2` matching the existing native valley-area scale. It is not a
   literature sediment-yield law. Mouth influence extends at most 512 m with
   `(1-d/512)^2`, taking the maximum rather than summing multiple samples of
   the same river. A connecting segment must remain within 64 m of the shore
   and cannot cross a lake. These distance/kernel choices are explicit
   finite-support closures, not measured longshore transport. No conservation
   claim is made for this index. Actual sediment flux and sediment inventory
   are necessary for an evolved beach or delta.

## Classification in evaluation order

Marine posts are `none` with missing-bed flag, never treated as intertidal
bed at z=0. Lakes are `none`. Only land within 512 m of the native shoreline
is considered. Wave forcing comes from the nearest reachable marine-touching ground anchor.
Slope is now computed at the actual classified post, not borrowed from that
anchor. All material classes require the width and elevation envelope
(`-TR/2<=z<=TR/2+R2`), without a first-post exception. A steep rock-face
annotation may extend above this envelope only at a face-adjacent marine post.

Native channel beds (codes 1/5), lakes and marine posts are not coastal
substrate. Ground must connect by four-neighbour non-water/non-channel posts
above low water to some ground above high water. Detached low bars with no
supratidal support are excluded with bit512. Small islands with a supratidal
core remain valid. This is a shore-attached product restriction, not a law
that real detached tidal shoals cannot exist. No native flags are relabelled.
A deterministic land-path search then links each candidate to a marine-touching
anchor within 512 m. It cannot jump across water as Euclidean nearest-anchor
selection could. Ground below high water must additionally be reachable from
marine contact through ground no higher than high water. This static access
check excludes low depressions behind supratidal barriers; it does not solve
tides, change water-body identity or use a marine surface as bed depth.

Steep => cliff. Else absent mobile-material evidence => rocky platform.
Mobile evidence is exposed shale/sandstone/conglomerate/coal-bearing cover or
unconsolidated material, or `S>=0.5`. This is inferred sediment availability,
not a measured sediment layer. With mobile evidence and effective mouth area
>=1 km2, shore distance <=32 m and slope <=3 degrees => delta-front potential.
Else `RTR>15`, local slope <=3 degrees and elevation within the tidal frame
=> tidal flat, or salt-marsh potential above the Balke et al. (2016) pioneer
boundary and no higher than MHW. There is no tidal-flat label above MHW. Marsh additionally carries a missing-ecology flag. No plant species,
salinity, colonization time or biomass is fabricated. Remaining mobile shores
are gravel-compatible at >=10 degrees or conglomerate exposure, otherwise
sand-compatible. Every mobile class is flagged as inferred substrate.

Dalrymple, Zaitlin and Boyd (1992) provide a wave/tide/river conceptual framework,
not the missing numerical thresholds for this task's eight-way inverse map.
No claim is made that the area proxy or delta eligibility threshold comes
from that paper. All three forcings are exposed honestly; a calibrated river
versus wave sediment budget is future work.

## Output contract

Each `w_RX_RY/coast.json` carries chart, seed, config, upstream window and the
complete context digest list. Lattice n=545, origin=(512*rx-16,512*ry-16),
32 m posts, row-major with x east/y north. All binary fields are little endian.

| Field | Type | Meaning |
|---|---|---|
| shore_class | u8 | 0 none, 1 cliff, 2 platform, 3 gravel, 4 sand, 5 flat, 6 marsh potential, 7 delta potential |
| coast_width_m | f32 | horizontal envelope; zero for none/cliff; censor flag at 512 m |
| sediment_supply_index | f32 | potential index 0..1; not mass or volume |
| wave_height_m | f32 | directional local wind-sea Hs proxy |
| shore_slope | f32 | local tangent at the actual classified post |
| coast_flags | u16 | bitmask below |

Bits: 1 inferred substrate; 2 marsh potential; 4 delta potential; 8 width
censored; 16 coarse T0 fetch used; 32 fetch capped; 64 marine bed missing;
128 assumed tide; 256 climate wind treated as U10; 512 unsupported shore
connection; 1024 channel bed; 2048 ground below high water but tidally isolated. Diagnostics are zero away
from the evaluated coast band. A zero class is not a missing-data sentinel;
read flags and the artifact's limitations. Class shares use half-open 512x512
cores so shared edges are not double counted.

## Unreal hookup proposal, not implemented

`t3_prep` would take an optional `--coast-dir` context. Verify chart, spacing,
seed, native window digest and common context before reading these six fields.
Read numeric fields at the existing terrain sample coordinates. Sample class
and bitfields categorically, never interpolate codes. Use sand/gravel/rock
material variants and the width field only as a bounded material mask on the
existing terrain and signed shoreline. Keep unknown/inferred status in the
asset receipt. Do not paint a marsh plant population from potential code 6
without the ecology lane, lower sea-plane terrain from a surface-only field,
or modify collision heights. Any future generated geometry must use a native
sediment/thickness or bathymetry model. No client or t3 file changed here.

Integration estimate: material-only adapter 0.5-1 day plus one native Unreal
inspection; physical intertidal geometry and sediment stocks 3-5 days after
river-sediment and marine-bed inputs exist. These estimates are proposals,
not authorization for additional work.


## Arc 2 marsh-source correction and remaining limits

Balke et al. (2016), Figure 5, give `Pio_h(cm) = -108.23 log10(MTR(cm)) +
163.21`, relative to MHW. For an assumed mean tide range of 2 m symmetric
about the sea datum, MHW=1 m and the pioneer lower bound is 0.141695 m.
The fit has R²=0.39 and is an empirical screening relation, not a species or
colonization model. [Primary paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015WR018318).
It does not justify an arbitrary higher cutoff to lower site5's class share.
The configured 2 m effective tide is used as MTR here and as a spring-range
proxy in RTR; the record has neither tidal constituent series nor separate
mean/spring datums. This approximation remains explicit.

RTR>15 screens wave-versus-tide regime; it alone does not establish marsh
substrate, salinity, species or recruitment. Wind-only exposure omits swell
and storms. Local low slopes and the pioneer elevation band still leave many
native constructed fan/bar posts eligible. A large area share can therefore
remain, or increase, after correcting geometry. Shares over classified land
area must not be presented as percentages of the volcanic coastline length.

The reviewer-identified rings in site5 are predominantly fixed=0, with native
over=2 fan and over=3 beach, on constructed surfaces from 0 to 1.8 m above
sea level. Some detached bars never exceed high tide and are rejected in v2.
Others contain supratidal ground, so their own intertidal shores remain.
Neither removing a real input island nor asserting its vegetation is justified.
Inspect the water-coloured same-crop figures and source-bound component audit
before integrating. Native bars and the missing marsh ecological inputs remain
upstream scope, not a reason to tune the classifier to a target percentage.
