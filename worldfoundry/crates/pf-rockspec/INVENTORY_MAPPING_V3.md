# Inventory to pf-rockspec v3 mapping

Canonical mineral IDs and dry-solid mass fractions are pf-rockspec IDs. Snapshot source SHA256: `241aaa6336a888d507d0ff80b0d806bbab60cbe8cd979d4f329f4fa775918a3c` (`pf-inventory/src/petrology.rs`). The source was read only. The verification script checks its density-match vocabulary against this mapping.

`w_i = phi_i rho_i / sum(phi_j rho_j)`. Supply and record every producer density in kg/m3. `phi` is solid volume fraction, not bulk volume including pores. A supplied species split is a dry-mass split within that group. No producer density is replaced silently.

| Producer ID (solid volume) | Supported pf-rockspec ID(s) | Required choice |
|---|---|---|
| silicic_glass | felsic_glass_low_absorption, volcanic_glass | Explicit body-admitted species/optical population split |
| quartz | quartz | Identity mapping; source density receipt |
| k_feldspar | k_feldspar | Identity mapping; source density receipt |
| plagioclase | plagioclase | Identity mapping; source density receipt |
| biotite | biotite | Identity mapping; source density receipt |
| clinopyroxene | augite | Explicit body-admitted species/optical population split |
| olivine | olivine | Identity mapping; source density receipt |
| amphibole | hornblende | Explicit body-admitted species/optical population split |
| calcite | calcite | Identity mapping; source density receipt |
| dolomite | dolomite | Identity mapping; source density receipt |
| clay | illite, smectite, kaolinite | Explicit body-admitted species/optical population split |
| muscovite | muscovite | Identity mapping; source density receipt |
| halite | halite | Identity mapping; source density receipt |
| gypsum | gypsum | Identity mapping; source density receipt |

Group candidates are the supported catalogue subset, not a claim of mineralogical synonymy. Clinopyroxene must be explicitly qualified as augite; amphibole as hornblende. Clay must resolve to the admitted illite/smectite/kaolinite fractions. Silicic glass must resolve its absorption population (pale felsic glass versus the existing dark glass proxy). Missing evidence is an error, not an average.

The next arc joins the actual inventory modes and provenance to body realizations. This mapping does not overwrite those modes or change pf-inventory. pf-fracture consumes the same canonical mass modes, mechanics in Pa and Pa sqrt(m), lengths in m, and friction angles in radians.
