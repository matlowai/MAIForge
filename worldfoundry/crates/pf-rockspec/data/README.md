All catalogue estimates are project priors. Bibliography is embedded in
catalog_v1.json. No numeric endpoint is represented as a transcribed measurement.
Units distinguish linear reflectance, dimensionless quantities, Mohs rank,
metres, radians, degrees, kg/m3, Pa and Pa sqrt(m).

Diagnostic and swatch inputs are frozen manufactured states explicitly allowed
by this brief. Their bare numbers are exact diagnostic selections, not measured
natural properties. Coating coverage, thickness and elapsed history are prescribed
fixtures, never output from weathering-rate simulation. Temperature 293.15 K,
water film 100 micrometres, layer thickness 100 micrometres and coat grain 10
micrometres are such selections. Unknown age is null. A coating mineral comes
from the weathering-family products but its quantity/history is not inferred.

The 128 cases include all 43 rock keys in dry and wet states, 34 weathered states,
and 8 targeted states: parallel/cross cleavage, abrasion, altered buried bulk on
a fresh cut, micro-footprint glass, zero-thickness layer, partial wetness and
rotated fabric. All three requested footprints occur. The 387-cell sheet provides
the full Cartesian product separately.

Glass hydration products, grus, dissolved species, vermiculite, leucoxene,
iddingsite and clay residue are pathway descriptors, not undeclared mineral-mode
IDs. If a consumer wishes to realize one as composition it must supply a known
mineral assemblage or extend the catalogue. Unknown pathways never default to
quartz, basalt or air.

## Version 2

`catalog_v2.json` adds declared connected-domain and structure priors.
`diagnostic_inputs_v2.json` and `swatch_inputs_v2.json` retain 128 and 387
versioned cases. V1 files stay immutable. See ../MODEL_V2.md for generation,
parameter units, exposure fixtures, source qualification and limits.
