//! Lane SHARP (2026-10-16, opt-in `ct_sharp=1`): 2 m re-profile of T1b's cliff-and-talus staircase (pf-region lane STEP /
//! STEP2, rule ct_remap_v1 + ct_remap_local_v1). LABEL: pre-signoff, fidelity tier Preview, derived, never canon; a
//! construction (no SimTime, no keyed draws).
//!
//! Why: T1b remaps bare-rock posts with a monotone staircase F_c of the pre-staircase surface h (published as ct_base_mm):
//! per period [s0, s1] the lower part rises at the talus tangent tb, the upper part at the cliff tangent ta. At 32 m the cliff
//! band is narrower than one post, so the B1 bicubic of the published heights smears each riser over 32-64 m (50-70 deg, not
//! the class angle). Here the same F_c is applied to the 2 m surface instead.
//!
//! What is read (all published by T1b; nothing private to pf-region is re-derived): ct_base_mm (h), height_mm (z = F(h)),
//! ct_regime, ct_cliff_permille (c), exposed_litho (class), the header's ct class table, step period and the rivernet rockmass
//! bounds. Per bare-rock T1b post k (regime 1 or 2, c < 1):
//!   G_k = |grad h| by central differences over 32 m (one-sided at the lattice edge), as ct_remap_local_v1 computes it;
//!   ta_k from c/ta + (1 - c)/tb = 1/G (the per-post rock-mass multiplier of the cliff tangent is not published); where the
//!   implied multiplier falls outside the published rockmass bounds, ta is clamped into them and G re-derived from c;
//!   r = tb/G, q = ta/G (the slopes of F on the talus and cliff limbs);
//!   the post's own period boundary from its own limb: z - h = (m - 1)(h - s), m = r (talus: s = s0) or q (cliff: s = s1).
//! The ladder of post k: its own boundary (exact) plus the boundaries of the bare posts within `radius_posts`, merged where
//! closer than `merge_tol_m` (median), levels between the own boundary and h dropped, and the next level kept only if it is
//! consistent with the own limb (a talus post below its period's talus top, a cliff post above it); missing levels continue
//! at the regular period P. F_k is then T1b's remap function on that ladder: F_k(h_k) = z_k exactly.
//! At a 2 m post with pipeline height z2 (built on the B1 bicubic of the published staircase): h2 = z2 - S(x, y), S the
//! bicubic of the T1b shift field (z - h; 0 off bare rock), so the steps are not applied twice; the new height is
//! h2 + sum_k w_k (F_k(h2) - h2) over the four T1b posts of the cell (bilinear weights; F = identity on soil / mantled posts).
//! Continuous in position and height; equal at a shared tile-edge post from either tile (pointwise); T1b posts reproduce T1b's
//! height up to the 2 m pipeline's own change there.

use serde::Serialize;

pub const SHARP_RULE: &str = "ct_sharp_v1 (lane SHARP, pf-detail opt-in ct_sharp=1): T1b's cliff-and-talus remap F_c (ct_remap_v1 / ct_remap_local_v1) re-applied at 2 m. Per bare-rock T1b post: G from central differences of ct_base_mm, cliff tangent from c/ta + (1-c)/tb = 1/G (clamped to the published rockmass bounds x the class cliff tangent), the post's own period boundary from its own limb (talus: s0, cliff: s1), a ladder of the neighbours' boundaries (radius_posts, merged within merge_tol_m), the regular step period beyond; at a 2 m post the pipeline height minus the bicubic of the T1b shift (height_mm - ct_base_mm) is the pre-staircase height h2 and the new height is h2 + sum_k w_k (F_k(h2) - h2) over the cell's four T1b posts (bilinear); channel interiors and sea posts untouched; applied before the B9 fill and the B8 mean re-correction (F varies in space, so it can open pits; the tile mean is re-centred on the B1 base as without the key)";
pub const SHARP_NOTE: &str = "lane SHARP (opt-in ct_sharp=1): bare-rock posts re-profiled at 2 m with T1b's own staircase remap (rule in the recipe params ct_sharp.rule): near-vertical risers at the post's cliff tangent and benches/talus at the talus tangent instead of the B1 bicubic's smeared steps; soil posts away from bare rock unchanged; reads only published T1b ct fields";
pub const SHARP_MM_NOTE: &str = "i32 mm lane SHARP ct_sharp=1: change of the height by the 2 m staircase re-profile (before the B9 fill and the B8 mean re-correction); 0 where untouched";
pub const SHARP_PRE_RULE: &str = "ct_sharp_pre_v1 (lane TREAD, opt-in ct_sharp_stage=pre): the same 2 m re-profile (ct_sharp_v1) applied to each tile domain's initial surface (B1 cut + B4 closure + lobe relief; channel interiors and sea untouched, corridor banks re-profiled but fixed as before) BEFORE the B5 relaxation and the DE solve instead of after quantization. Riser mask = posts on the cliff limb (2 m regime 2). B5: riser posts are held (not crept or incised). DE: riser posts are held in the layered equilibrium (Picard) and the guard; in the threshold-failure step they are free with their own rock-mass cliff tangent (bilinear over the cell's bare T1b posts) as the failure threshold and no incoming-slope exemption (rock steeper than its cliff angle, e.g. 2 m detail stretched by the remap, falls); debris crossing a riser avalanches to its foot and rests at the talus repose slope (de_repose). Bench posts evolve under the unchanged DE rules. The post-stage step is not run; ct_regime / sharp_mm are the tile owner's domain values at the initial surface";
pub const SHARP_TREAD_RULE: &str = "ct_sharp_tread_v1 (lane SHARP2; on whenever ct_sharp=1 reads a T1b window carrying pf-region ct_tread's ct_tread_permille / ct_tread_tan_permille, else ct_sharp_v1 unchanged): T1b posts with a 3-limb profile (ft > 0) re-profile at 2 m with pf-region's ct_tread_v1 function: per period [s0, s1] (height T) climbing up, the lowest ft T is TREAD at tt (the bench on the top of the cliff below), the next fb T = (1 - ft - fc) T TALUS at the class tangent tb, the upper fc T CLIFF at ta. Per post: ft, tt and fc (ct_cliff_permille) are published, G from central differences of ct_base_mm as ct_sharp_v1, ta from ft/tt + fb/tb + fc/ta = 1/G (clamped to the published rockmass bounds x the class cliff tangent, G then re-derived), F-slopes p = tt/G (tread), r = tb/G (talus), q = ta/G (cliff). Own period boundary: a tread post gives s0 exactly (z - h = (p - 1)(h - s0)), a cliff post s1 exactly (as ct_sharp_v1); a talus post sits between two unknown boundaries (z - h = (r - 1)(h - s0) - T ft (r/p - 1)): its period is solved from the nearest neighbour boundary (exact tread / cliff estimates only, same radius and merge tolerance) above h taken as s1, else below h taken as s0, the first giving a consistent talus position; else T = the step period, else the period centring the post on its talus limb (counted). Ladder consistency (tread post: h on the tread of [own, up]; cliff post: h on the cliff of [dn, own]) and the regular-period extension as ct_sharp_v1 with the 3-limb breakpoints. T1b posts without a tread (2-limb fallback, massive classes, apron period) keep ct_sharp_v1. The 2 m application (h2 = z2 - bicubic shift, bilinear over the cell's four T1b posts), the post / pre stage and the riser mask (regime 2) are unchanged; tread posts are 2 m regime 3 (not held in the pre stage)";
pub const SHARP_TREAD_NOTE: &str = "lane SHARP2 (ct_sharp=1 on a window with pf-region ct_tread fields): bare-rock T1b posts with a 3-limb profile re-profiled at 2 m as tread (structural bench at the published tread tangent on the top of the cliff below), talus at the class repose and cliff at the post's cliff tangent (rule in block.json ct_sharp.window_build.tread.rule); 2 m ct_regime 3 = tread";
pub const SHARP_REGIME_TREAD_NOTE: &str = "u8 lane SHARP ct_sharp=1 + lane SHARP2 tread fields: 2 m regime from the re-profile: 0 soil/mantle (or less than half the bilinear weight on bare-rock T1b posts), 1 talus/bench at repose, 2 cliff/riser, 3 tread (structural bench at the published tread tangent) (the limb of F at the post's pre-staircase height, for the T1b post of the cell with the largest bare-rock weight)";
pub const SHARP_BLEND_PHASE_RULE: &str = "ct_sharp_blend_phase_v1 (lane PHASE, opt-in ct_sharp_blend=phase): at a 2 m post each bare-rock T1b corner of the cell re-profiles with its F on a BLENDED ladder instead of its own: every level L of the corner's ladder (extended by the regular step period) becomes the weighted mean of L (bilinear weight w_k) and each other bare corner's nearest level n_j (weight w_j x agree(|n_j - L|), agree = 1 within merge_tol_m, fading linearly to 0 at 2 merge_tol_m); the period [s0*, s1*] holding the pre-staircase height h2 is searched on the blended ladder; the knot fractions of F (u and (z - s0)/T at the tread -> talus and talus -> cliff breaks; 2-limb corners: no tread) are blended with the same weights (agreement at both boundaries); F is piecewise linear through the knots (s0*, s0*), (s0* + T u1, s0* + T z1), (s0* + T u2, s0* + T z2), (s1*, s1*). Corners whose ladders agree (the same levels a few m apart) thus share ONE F at the post (a bench is one flat limb across a mixed cell and the limb boundaries are one line) instead of a blend of four F values whose limbs straddle each other; corners that disagree by more than twice the tolerance keep their own F (on a cliff limb a boundary offset dh moves the surface by (q - 1) dh, so far-apart levels are not blended). The new height is h2 + sum_k w_k (Fk*(h2) - h2) over the bare corners (soil corners identity, as before); the 2 m regime is the limb of the heaviest bare corner's Fk* (0 when the bare weight < 0.5); the cliff tangent stays the bilinear mean. Continuous in position (no reference corner); F keeps its shape rules (identity at the blended boundaries, monotone, limb slopes within the corners' range); a post with one bare corner of positive weight is that corner's F exactly (T1b posts reproduced); a period walk that does not converge (6 moves or a flip) uses the corner's own F";
pub const SHARP_BLEND_NEAREST_RULE: &str = "ct_sharp_blend_nearest_v1 (lane PHASE diagnostic, opt-in ct_sharp_blend=nearest): at a 2 m post only the dominant bare corner's F (largest bilinear weight) re-profiles the post, h2 + wb (F_ref(h2) - h2); discontinuous at the T1b cell mid-lines; not a production surface";
pub const SHARP_BLEND_PERIOD_RULE: &str = "ct_sharp_blend_period_v1 (lane SPLIT-2M, opt-in ct_sharp_blend_period=1, needs ct_sharp_blend=phase): in the phase blend another bare corner j lends its knot fractions to corner k's period [s0, s1] only when j's nearest levels to s0 and s1 bound ONE period of j's own ladder (j's period holding (s0 + s1) / 2 is [near_j(s0), near_j(s1)]), not merely when j has some level near each end. Measured cause (lane SPLIT-2M, preview-client south wall on the vs=1 split window): a post A whose period is [-1347, 322] m (no riser bed in between at A, the bed there is under the 48 m riser threshold) and its neighbour B with [s_dn, s0, s1] = [-1347, -109, 322] m (the same bed IS a riser at B) agree at both of A's ends (B's s_dn and s1), so A's 1669 m period took B's knot fractions, which describe B's 431 m period: the knots moved by up to 90 m and the 2 m surface between A and B was lowered 20-29 m in a band along every such ladder jump (the bead rows). Fractions are only comparable on the same period. The level blend, the period walk and every other rule are unchanged; F is the identity at the blended period ends whatever the knots, so the result stays continuous in position and height (the condition is fixed per corner pair and period).";
pub const SHARP_SOIL_EDGE_RULE: &str = "ct_sharp_soil_edge_v1 (lane TREADSOIL, opt-in ct_sharp_soil_edge=1): in a T1b cell with a soil-mantled corner (T1b regime 0) the 2 m staircase no longer fades linearly with the bilinear soil weight ws across the whole cell (h2 + sum_k w_k (Fk - h2): a bench tilts toward the soil corner by its own cut depth (1 - p)(h2 - s0) over 32 m). A soil mantle on a bench is thin colluvium over the same flat bench, so on the tread limb the bare corners' staircase is kept whole (renormalised over the bare weight, h2 + sum_k w_k (Fk - h2) / (1 - ws)) up to the mantle edge ws = e0 = 0.5 (the 2 m regime boundary, which is unchanged) and the colluvium wedge h2 - F thickens from 0 there to the full T1b soil surface at the soil corner over the soil side of the cell (factor 1 / (1 - e0) on the w_k). The wedge stays a mantle only while its extra gradient |Fk - h2| / ((1 - e0) 32 m) stays within TW = 1.0 (sub-key ct_sharp_soil_edge_tw, 0 = no limit): e0 = min(0.5, 1 - |Fk - h2| / (32 m TW)) per corner and post (the whole edge up to a 16 m cut), so a soil corner far above or below the extended bench (|Fk - h2| >= 32 m) keeps the linear fade (without the limit a too-steep wedge relaxes in B5 and pulls soil T1b posts off their height: crop soil-post max error 9.7 -> 19.9 m). Per bare corner k the factor is min(1 / (1 - ws), 1 / (1 - e0 e_k)), e_k = 1 on its tread limb, falling linearly to 0 across its talus limb (continuity at the tread -> talus break), 0 on its cliff limb and for 2-limb corners (the cliff and talus fade as before; F - h2 = 0 at the period boundaries, so the jump of e_k there is continuous in height). ct_sharp_blend=phase: e_k from the corner's blended knots. Soil posts on or next to (8-neighbour) a fixed T1b post (channel corridor, sea: valley-floor and shore soil, not bench colluvium) are not soil corners for this rule (linear fade), so an extended bench never undercuts a channel corridor (F1). Cells without a soil corner, all T1b posts (soil: every w_k = 0; bare: ws = 0) and non-soil identity corners (steeper than cliff, no class, ladder failed: linear fade as before) are unchanged; the result is continuous (ws and the corners' F are edge-local on a cell edge)";
pub const SHARP_CARRY_STRIKE_RULE: &str = "ct_sharp_carry_strike_v1 (lane CARRY, opt-in ct_sharp_carry_strike=1, needs ct_sharp_stage=pre): the held riser posts are carried after B5 and after the DE equilibrium by a change that is consistent along the riser. Each riser post first takes the change of the least-squares plane through the non-riser posts of its 5 x 5 window (lane TREAD), then the carried changes are low-passed along the riser's contour: a Gaussian (sigma = ct_sharp_carry_strike_m, default 8 m; taps every 2 m to 2.5 sigma) over the carried posts met by tracing the contour both ways from the post (strike of the domain's pre-B5 surface smoothed by 3 m; bilinear over carried posts, a mask weight < 0.3 ends the band). Basis: a cliff band is one rock body; B5 / DE relaxation of the slopes around it translates the band, it does not shear it post by post. The plane fit's post-to-post scatter comes from each window holding a different set of talus / tread neighbours that moved metres (lane FLUTE: corr 0.985 with the riser B5 residual, 52 % of the riser along-strike high-pass), i.e. window sampling, not rock. The fall-line variation (talus below vs tread above moving differently) is kept: the path follows the contour. Riser failure in the DE failure step is unchanged. Deterministic, no draws";
pub const SHARP_LIMB_LP_RULE: &str = "ct_sharp_limb_lp_v1 (lane CARRY, opt-in ct_sharp_limb_lp=1, needs ct_sharp_stage=pre): a post whose staircase lands on the cliff limb (2 m regime 2 riser) is re-evaluated with its B4 closure replaced by the closure low-passed along the T1b contour (Gaussian sigma = ct_sharp_limb_lp_m, default 8 m, taps every 2 m to 2.5 sigma, path traced along the strike of the T1b bicubic surface, closure bilinear on the 2 m lattice). Basis: the closure is unresolved meso relief of the slope (keyed draws); the cliff limb's F slope (>> 1) multiplies it into post-to-post steps of the rock face (lane FLUTE: riser sharp_mm high-pass on closure high-pass corr 0.46, gain 2.25) that are not rock structure. A pure function of the post position (the closure field and the T1b surface are global), so shared posts agree across domains. Non-riser posts keep their own closure. Deterministic, no new draws";
/// Lane CARRY: default along-contour sigma (m) of ct_sharp_carry_strike and ct_sharp_limb_lp.
pub const CARRY_STRIKE_M: f64 = 8.0;
pub const SHARP_AA_RULE: &str = "ct_sharp_aa_v1 (lane CARRY2, opt-in ct_sharp_aa=1, needs ct_sharp_stage=pre): a 2 m post's staircase height is the MEAN of the staircase surface over the post's footprint (a square of ct_sharp_aa_w posts, default 1 = the 2 m cell, sampled 4 x 4), not its value at the post point. The staircase surface at a footprint point is the post's F-map applied there (SharpField::apply_ta at the point, its own T1b cell weights and bicubic shift) to the pre-staircase surface there (bilinear over the pre-staircase posts, a pure function of position: B1 cut + B4 closure + lobe relief, closure 0 on corridor / sea posts; shifted by the post's own closure swap of ct_sharp_limb_lp). Written as zn = mean_p F(z_p) - sigma (mean_p z_p - z_c) with sigma = dF/dz at the post (central difference, 1 cm): where F is linear over the footprint the post keeps its point value exactly (to the positional curvature of the blend), so only posts whose footprint holds a knot line (cliff foot, cliff top, tread edge, period boundary) change, by at most the footprint mean of the kink. The 2 m regime, cliff tangent and riser mask stay the point values; fixed posts (corridor banks) keep the point value. Basis: a 2 m DEM post stands for its footprint; the staircase's knot lines are sharp (a 50 m riser 3 posts wide at s2_11), and point-sampling them on the lattice aliases every knot line that runs oblique to the lattice into post-to-post steps (lane CARRY2: along-strike high-pass of the initial staircase surface on s2_11 lowest-riser posts 0.17 m where the contour runs within 5 deg of a lattice axis vs 0.60-0.67 m at 15-45 deg; the staircase change carries 97 % of it, its input h2 3 %), which the renderer's triangles show as vertical hatching. The footprint mean puts each knot line's sub-post position into the posts (anti-aliasing), so the bilinear / triangle surface between posts follows a straight knot line. Deterministic, no draws";
pub const SHARP_RISER_NOFLUX_RULE: &str = "ct_sharp_riser_noflux_v1 (lane CARRY2, opt-in ct_sharp_riser_noflux=1, needs ct_sharp_stage=pre): in the B5 relaxation the held riser posts (2 m regime 2 on the cliff limb, not corridor / sea) are a ZERO-FLUX boundary of the nonlinear creep: a 4-neighbour creep edge with a riser post at either end carries no flux (before: the held riser kept its height while every creep step exchanged flux with its free neighbours, a fixed-height source below the cliff foot and sink above the crest, the flux at the critical-slope cap fmax x the 5-15 m height step per post, applied to the free post only). Stream power, drainage area and the riser hold itself are unchanged; the held-riser carry follows as before. Basis: the riser is bare rock, it supplies no creeping regolith and is not a reservoir; mass enters the foot from the cliff only through the DE failure step (rockfall to talus at repose). Lattice: a riser oblique to the 2 m lattice has a 1-post jog in its mask every 1 / tan(angle) posts; with the fixed-height exchange each free post took flux in proportion to its number of riser 4-neighbours, so every jog grew a knob at the foot and a notch at the crest (lane CARRY2: s2_11 lowest riser, periodic knobs every 4-6 posts appear in B5, the staircase surface before B5 is smooth along the riser). Deterministic, no draws";
pub const SHARP_CARRY_FRINGE_RULE: &str = "ct_sharp_carry_fringe_v1 (lane CARRY2, opt-in ct_sharp_carry_fringe=1, needs ct_sharp_carry_strike=1): after B5 the riser FRINGE (free posts, not riser, not corridor / sea / ring, within FRINGE_POSTS = 2 posts, Chebyshev, of a held riser post: the cliff foot apron and the crest) takes its B5 change blended toward its low-pass along the contour with weight (3 - d) / 2 at Chebyshev distance d (1 next to the riser, 0.5 at 2 posts, 0 beyond: no step at the fringe edge), the low-pass as the held-riser carry (ct_sharp_carry_strike: Gaussian sigma ct_sharp_carry_strike_m, path along the strike of the pre-B5 surface smoothed 3 m, bilinear over fringe posts, mask weight < 0.3 ends the band), before the riser carry reads its 5 x 5 plane-fit sources. Basis: B5's creep exchanges flux with each held riser post across 4-neighbour edges, so on a riser oblique to the 2 m lattice (a 1-post jog in its mask every 1 / tan(angle) posts) each foot / crest post takes 0, 1 or 2 riser edges' worth of supply depending on where the jog falls: a knob at every jog of the foot and a notch at the crest (lane CARRY2: s2_11 lowest riser, knobs every 4-6 posts appear in B5; the staircase before B5 is smooth along the riser). The supply per unit length of one cliff band is a property of the band, not of the lattice phase; the low-pass keeps the band's mean exchange (riser height, foot apron volume) and its fall-line profile (the path follows the contour). Deterministic, no draws";
/// Lane CARRY2: riser fringe width (posts) of ct_sharp_carry_fringe.
pub const FRINGE_POSTS: i64 = 2;
/// Lane CARRY2: default footprint width (posts) of ct_sharp_aa and the sub-samples per axis.
pub const AA_W_POSTS: f64 = 1.0;
pub const AA_SAMPLES: usize = 4;
pub const LIMB_LP_M: f64 = 8.0;
pub const SHARP_TALUS_ROCK_RULE: &str = "ct_sharp_talus_rock_v1 (lane STEEP-TALUS, opt-in ct_sharp_talus_rock=1, needs ct_sharp_stage=pre): the 2 m regime follows the surface it labels. A 2 m post takes the limb (tread / talus / cliff) of its heaviest bare T1b corner's F at its own pre-staircase height, but the surface between reproduced T1b posts climbs whatever those posts climb: where neighbouring T1b talus posts sit on different periods of the staircase (their cliff falls between the 32 m posts and their ladders disagree beyond the phase-blend tolerance) the 2 m surface is a 50-60+ deg ramp whose every post is on a corner's talus limb (lane STEEP-TALUS: region R 2 m regime-1 posts > 50 deg at 8 m sit 82 % in T1b cells whose bilinear T1b height already climbs > 50 deg, gentle talus 4 %). Coarse debris cannot stand there (repose 34-40 deg; the DE talus repose 0.75): the ramp is the T1b cliff band smeared over the cell, i.e. rock. Rule, on each tile domain's pre-staircase surface (B1 cut + B4 closure + lobe relief + the 2 m staircase, before B5 / DE): a regime-1 post whose surface is steeper than tan(theta) over 8 m both up- and downslope (per axis the minmod of the one-sided 4-post (8 m) differences; posts closer than 4 posts to the domain ring are not tested; theta = 45 deg, sub-key ct_sharp_talus_rock_deg) is relabelled 2 (rock). It does not join the riser mask: B5 and the DE treat it as the bare rock it is (DE's slope exposure in the layered equilibrium, failure at the span critical slope with the incoming-slope exemption, as any bare rock); it is not a de_talus_mantle / de_talus_runoff talus post. The minmod keeps the talus fringe at a cliff foot / crest (steep on one side only) talus. Heights are not changed by the relabel itself (only through the B5 / DE treatment of the relabelled posts). Regime-0 (soil) and regime-3 (tread) posts are not relabelled";
pub const SHARP_ENDS_RULE: &str = "ct_sharp_ends_v1 (lane PERIOD-ENDS-2, opt-in ct_sharp_ends=1, needs a T1b window with pf-region ct_ends's fields): each bare-rock T1b post with a published period (ct_s0_mm / ct_s1_mm, kind 1-6) whose height is F(ct_base_mm) on its published knots (kind without the +128 flag) re-profiles at 2 m with its EXACT F instead of the one ct_sharp_v1 / ct_sharp_tread_v1 infer: knot fractions u = (h - s0) / T and f = (F - s0) / T at the tread -> talus (ct_knot_h1_mm / ct_knot_z1_mm; none when h1 = s0) and talus -> cliff (ct_knot_h2_mm / ct_knot_z2_mm) breaks, F piecewise linear (s0, s0) -> knot 1 -> knot 2 -> (s1, s1), the same fractions in every period of the ladder; ladder = the post's own [s_dn, s0, s1, s_up] (ct_s_dn_mm / ct_s_up_mm where published), beyond its ends the regular step period (counted as s_dn / s_up missing); apron posts (kind 5, s0 = the hillslope base): [s0, s1] at and above s0, [s_dn, s1] below (pf-region's remap below the base), same fractions; massive posts (kind 6): identity outside [s0, s1]. The cliff tangent (DE riser failure threshold) is S x (1 - f2) / (1 - u2) with S = ct_grad_permille (the remap's own gradient; no central differences, no clamp). No own-boundary estimate from the shift, no talus period solve, no neighbour ladder merge for these posts. Posts flagged +128 (held by the tread base-level / drainage rule or moved after the remap: bank, fill, coast) and bare posts without a published period keep the inferred corner of ct_sharp_v1 / ct_sharp_tread_v1 (built from the same neighbour estimates as without the key). The 2 m application (h2, phase blend, soil edge, pre stage, riser mask) is unchanged";
pub const SHARP_ENDS_LABEL_RULE: &str = "ct_sharp_ends_label_v1 (lane PERIOD-ENDS-2 diagnostic, opt-in ct_sharp_ends_label=1, needs ct_sharp_ends=1): the 2 m regime is the limb of the BLENDED F at the post's pre-staircase height instead of the heaviest bare corner's limb: blended F-slope sb = 1 + sum_k g_k w_k (F'_k(h2) - 1) / sum_k w_k over the cell's bare corners (F'_k the slope of corner k's (phase: blended) F segment at h2, g_k the soil-edge factor), labelled cliff (2) where sb >= sqrt(r q), tread (3) where the heaviest corner has a tread and sb <= sqrt(p r), else talus (1), with p, r, q the heaviest corner's tread / talus / cliff F-slopes (its limb where it has no talus or cliff limb); 0 where the bare weight < 0.5 as before. A hidden riser between T1b posts on different periods (the lower corner on its cliff limb, the upper on its talus) is then labelled by the slope the blend actually builds";
pub const SHARP_CELL_RISER_RULE: &str = "ct_sharp_cell_riser_v1 (lane CELL-RISER, opt-in ct_sharp_cell_riser=1, needs ct_sharp_blend=phase and ct_sharp_ends=1): a cliff band is one outcrop of the hard unit; its foot and top follow the bed, they do not fan out across a 32 m cell. pf-region's remap gives every T1b post its own talus -> cliff break (knot fractions from the post's own gradient S), so neighbouring posts on the SAME period can put the cliff foot 50-80 m apart (lane PERIOD-ENDS-2: s2_09 knot f2 0.32-0.80 over one cell; 7-12 % of same-period neighbour pairs > 20 m); reproducing both posts exactly while blending their knots across the cell smears the climb between them into a 55-62 deg ramp labelled talus. Rule: a bare T1b post with an exact corner (ct_sharp_ends; not apron, not massive) whose cliff foot z2 = s0 + T f2 lies more than D (default 20 m, sub-key ct_sharp_cell_riser_m) above the lowest cliff foot of its same-period 8-neighbours (s0 and s1 both within merge_tol_m; exact, non-apron, non-massive) re-profiles with that neighbour's knot fractions (u1, f1, u2, f2 on its own ladder; the lowest foot = where the band is fully exposed, the higher feet are the same band buried by more talus) and keeps its height EXACTLY by a pseudo-base shift dh = F*^-1(z) - h (F* monotone on the post's ladder, so the inverse is unique; a post above the lowest foot thus sits on the cliff limb of F*). At a 2 m post the corners' F are evaluated at h2 + dh(x), dh(x) = the cell's four per-post shifts (0 at posts without one) blended with smoothstep weights 3t^2 - 2t^3 per axis: it passes through every post's dh, is C1 across the T1b cell edges and never overshoots (sub-key ct_sharp_cell_riser_dh: smooth default; bilinear, tried first, leaves a slope kink on every cell edge, 2 m curvature at the cell-edge phase +13-25 %; bicubic Catmull-Rom is C1 but overshoots between alternating shifts and rebuilds steep talus), and the height is h2 + sum_k g_k w_k (Fk*(h2 + dh) - h2) (phase blend, soil edge as before). The climb between a corner below the foot and one above it is then talus at the F* talus slope up to ONE break contour (the level h2 + dh = break of F*) and the cliff limb above it at the F* cliff slope, instead of a ramp: the shift moves where the break sits, not how steep the limbs are. T1b posts are reproduced (F*(h + dh) = z; knots in mm), every field is a per-T1b-post value interpolated continuously (knots by the phase blend's bilinear weights, dh as above; on a cell edge only its two posts count), F* is monotone per post; posts that move nothing (foot within D of the lowest, or the post is the lowest) have dh = 0 and their own knots. The cliff tangent (DE riser threshold) of a re-profiled post is the donor's S q";
pub const CELL_RISER2_RULE: &str = "ct_sharp_cell_riser2_v1 (lane CELL-RISER-2, opt-in sub-keys of ct_sharp_cell_riser): (a) ct_sharp_cell_riser_donor: min (default, v1) compares the post's cliff foot with the LOWEST same-period neighbour foot, so one noisy low neighbour re-profiles a whole cell with a large shift; median = the lower median of the same-period neighbour feet (the band's typical exposure around the post: a post is re-profiled only when its foot is an outlier above most of its neighbours, and takes the typical neighbour's knots, a smaller shift); pair = the lowest foot that a second same-period neighbour confirms (its foot within D of it), so a single low outlier is not a donor. (b) ct_sharp_cell_riser_cap=<cells>: the shift moves the post's break contour horizontally by |dh| / S (S = the remap gradient ct_grad_permille); the v1 premise is that the higher foot is the same band within one cell, so a post is re-profiled only when |dh| <= cap x 32 m x S (posts beyond keep their own knots, dh 0); 0.5 keeps the break in the post's own half of every cell edge. (c) ct_sharp_cell_riser_dh=monotone: the per-post shift dh at a 2 m post by a tensor-product monotone cubic Hermite (Fritsch-Carlson class, Fritsch-Butland harmonic-mean tangents: C1, monotone between posts, no overshoot of the two posts' values, tangents non-zero at posts where the data are monotone, unlike the smoothstep's zero slope at every post)";
pub const SHARP_REGIME_NOTE: &str = "u8 lane SHARP ct_sharp=1: 2 m regime from the re-profile: 0 soil/mantle (or less than half the bilinear weight on bare-rock T1b posts), 1 talus/bench, 2 cliff/riser (the limb of F at the post's pre-staircase height, for the T1b post of the cell with the largest bare-rock weight)";

#[derive(Clone, Debug, Serialize)]
pub struct SharpParams {
    pub rule: String,
    /// Boundary estimates of neighbouring posts closer than this (m) are one level.
    pub merge_tol_m: f64,
    /// Chebyshev radius (T1b posts) of the ladder neighbourhood.
    pub radius_posts: i64,
    /// Lane TREAD (opt-in `ct_sharp_stage=pre`): None = the post stage (after quantization; every byte as lane SHARP).
    /// Some(rule) = the pre stage (SHARP_PRE_RULE): the staircase is applied to each tile domain's initial surface before the B5
    /// relaxation and the DE solve, with the riser posts protected.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stage: Option<String>,
    /// Lane SHARP2 (opt-in `ct_sharp_consensus=1`, 3-limb posts only): a tread / cliff post's own period boundary is replaced by
    /// the median of its merge cluster (the neighbours' estimates within merge_tol_m), so neighbouring corners share the level
    /// and the bilinear blend does not tilt the treads; T1b posts are then reproduced only to (1 - p) x the cluster spread.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub t3_consensus: bool,
    /// Lane PHASE (opt-in `ct_sharp_blend=phase|nearest`): None = the bilinear blend of the four corners' F values (every byte
    /// as before). Some(SHARP_BLEND_PHASE_RULE) = one F per 2 m post from the corners' blended period knots; Some(SHARP_BLEND_NEAREST_RULE)
    /// = the dominant corner's F only (diagnostic, discontinuous at the cell mid-lines).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub blend: Option<String>,
    /// Lane TREADSOIL (opt-in `ct_sharp_soil_edge=1`): None = the linear soil fade over the cell (every byte as before);
    /// Some(SHARP_SOIL_EDGE_RULE) = benches next to a soil corner stay whole up to the mantle edge (ws = 0.5).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub soil_edge: Option<String>,
    /// Lane TREADSOIL (opt-in `ct_sharp_soil_edge_tw=<tangent>`, needs ct_sharp_soil_edge=1): the wedge's extra-gradient limit
    /// (None = SOIL_EDGE_TW 1.0; 0 = no limit, diagnostic).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub soil_edge_tw: Option<f64>,
    /// Lane STEEP-TALUS (opt-in `ct_sharp_talus_rock=1`, pre stage only): None = the 2 m regime is the heaviest corner's limb
    /// (every byte as before); Some(SHARP_TALUS_ROCK_RULE) = talus posts steeper than talus can stand at 8 m are relabelled rock.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub talus_rock: Option<String>,
    /// Lane STEEP-TALUS (opt-in `ct_sharp_talus_rock_deg=<deg>`, needs ct_sharp_talus_rock=1): the threshold angle (None = TALUS_ROCK_DEG).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub talus_rock_deg: Option<f64>,
    /// Lane PERIOD-ENDS-2 (opt-in `ct_sharp_ends=1`, needs a window with pf-region ct_ends's fields): None = the inferred ladders
    /// (every byte as before); Some(SHARP_ENDS_RULE) = each bare T1b post's F from its published period ends and knots.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ends: Option<String>,
    /// Lane PERIOD-ENDS-2 (opt-in `ct_sharp_ends_label=1`, needs ct_sharp_ends=1): the 2 m regime is the limb of the blended F
    /// (its slope at h2) instead of the heaviest corner's limb (SHARP_ENDS_LABEL_RULE).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ends_label: Option<String>,
    /// Lane CELL-RISER (opt-in `ct_sharp_cell_riser=1`, needs ct_sharp_blend=phase and ct_sharp_ends=1): None = every byte as
    /// before; Some(SHARP_CELL_RISER_RULE) = same-period posts whose cliff foot sits more than D above a neighbour's take that
    /// neighbour's knots and keep their height by a pseudo-base shift (one break contour per cell instead of a ramp).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_riser: Option<String>,
    /// Lane CELL-RISER (opt-in `ct_sharp_cell_riser_m=<m>`, needs ct_sharp_cell_riser=1): D (None = CELL_RISER_D_M).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_riser_m: Option<f64>,
    /// Lane CELL-RISER (opt-in `ct_sharp_cell_riser_dh=bilinear|bicubic|smooth`, needs ct_sharp_cell_riser=1): how the per-post
    /// shift dh is interpolated at a 2 m post (None = CELL_RISER_DH_DEFAULT).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_riser_dh: Option<String>,
    /// Lane CELL-RISER-2 (opt-in `ct_sharp_cell_riser_donor=median|pair`, needs ct_sharp_cell_riser=1): which same-period
    /// neighbour's foot the post is compared with and takes its knots from (None = the lowest foot, CELL_RISER2_RULE).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_riser_donor: Option<String>,
    /// Lane CELL-RISER-2 (opt-in `ct_sharp_cell_riser_cap=<cells>`, needs ct_sharp_cell_riser=1): a post is re-profiled only when
    /// its shift moves the break contour by at most this many T1b cells, |dh| <= cap x 32 m x S (None = no cap).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_riser_cap: Option<f64>,
    /// Lane CARRY (opt-in `ct_sharp_carry_strike=1`, pre stage only): None = each held riser post is carried by its own 5 x 5
    /// plane fit (every byte as before); Some(SHARP_CARRY_STRIKE_RULE) = the carried change is low-passed along the riser contour.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub carry_strike: Option<String>,
    /// Lane CARRY (opt-in `ct_sharp_carry_strike_m=<m>`, needs ct_sharp_carry_strike=1): the along-contour Gaussian sigma (None = CARRY_STRIKE_M).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub carry_strike_m: Option<f64>,
    /// Lane CARRY (opt-in `ct_sharp_limb_lp=1`, pre stage only): None = the staircase reads the post's own B4 closure (every byte
    /// as before); Some(SHARP_LIMB_LP_RULE) = cliff-limb posts read the closure low-passed along the T1b contour.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub limb_lp: Option<String>,
    /// Lane CARRY (opt-in `ct_sharp_limb_lp_m=<m>`, needs ct_sharp_limb_lp=1): the along-contour Gaussian sigma (None = LIMB_LP_M).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub limb_lp_m: Option<f64>,
    /// Lane CARRY2 (opt-in `ct_sharp_aa=1`, pre stage only): None = each staircase post is the staircase at the post point (every
    /// byte as before); Some(SHARP_AA_RULE) = the staircase's mean over the post's footprint (anti-aliased knot lines).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aa: Option<String>,
    /// Lane CARRY2 (opt-in `ct_sharp_aa_w=<posts>`, needs ct_sharp_aa=1): the footprint width in posts (None = AA_W_POSTS).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aa_w: Option<f64>,
    /// Lane CARRY2 (opt-in `ct_sharp_riser_noflux=1`, pre stage only): None = B5 creep exchanges flux with the held riser posts
    /// (fixed-height source / sink; every byte as before); Some(SHARP_RISER_NOFLUX_RULE) = no creep flux across a riser edge.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub riser_noflux: Option<String>,
    /// Lane CARRY2 (opt-in `ct_sharp_carry_fringe=1`, needs ct_sharp_carry_strike=1): None = the riser fringe keeps its own B5 change
    /// (every byte as before); Some(SHARP_CARRY_FRINGE_RULE) = the fringe's B5 change low-passed along the contour.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub carry_fringe: Option<String>,
    /// Lane SPLIT-2M (opt-in `ct_sharp_blend_period=1`, needs ct_sharp_blend=phase): None = every byte as before;
    /// Some(SHARP_BLEND_PERIOD_RULE) = a corner lends its knots in the phase blend only on its own matching period.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub blend_period: Option<String>,
}

impl SharpParams {
    /// True for the pre stage (`ct_sharp_stage=pre`).
    pub fn pre(&self) -> bool { self.stage.is_some() }
}

impl Default for SharpParams {
    fn default() -> Self { SharpParams { rule: SHARP_RULE.to_string(), merge_tol_m: 6.0, radius_posts: 3, stage: None, t3_consensus: false, blend: None, soil_edge: None, soil_edge_tw: None, talus_rock: None, talus_rock_deg: None, ends: None, ends_label: None, cell_riser: None, cell_riser_m: None, cell_riser_dh: None, cell_riser_donor: None, cell_riser_cap: None, carry_strike: None, carry_strike_m: None, limb_lp: None, limb_lp_m: None, aa: None, aa_w: None, riser_noflux: None, carry_fringe: None, blend_period: None } }
}

/// One class of the T1b ct table: (mode 0 none / 1 massive / 2 layered, talus tangent, cliff tangent).
#[derive(Clone, Copy, Debug)]
pub struct CtClass { pub mode: u32, pub tb: f64, pub ta: f64 }

/// The published T1b ct inputs over a window lattice (row-major n x n from (gx0, gy0)).
pub struct CtInputs<'a> {
    pub n: usize,
    pub gx0: i64,
    pub gy0: i64,
    /// ct_base_mm (m) and height_mm (m).
    pub h: &'a [f64],
    pub z: &'a [f64],
    pub regime: &'a [u8],
    pub cliff_permille: &'a [u16],
    pub litho: &'a [u8],
    pub classes: Vec<CtClass>,
    pub step_m: f64,
    pub rockmass: (f64, f64),
    /// Lane SHARP2: ct_tread_permille (ft) and ct_tread_tan_permille (tt) when the window carries pf-region ct_tread's fields.
    pub tread: Option<(&'a [u16], &'a [u16])>,
    /// Lane PERIOD-ENDS-2: pf-region ct_ends's fields (read only with ct_sharp_ends=1).
    pub ends: Option<&'a crate::input::CtEnds>,
}

/// Lane SHARP2: the 3-limb profile of one T1b post (pf-region ct_tread_v1): F-slope p = tt / G of the tread limb and the
/// height fractions (tread ft, talus fb, cliff fc) of the period; ft/p + fb/r + fc/q = 1.
#[derive(Clone, Copy, Debug)]
pub struct Tread3 { pub p: f64, pub ft: f64, pub fb: f64, pub fc: f64 }

impl Tread3 {
    /// The u = (h - s0)/T breakpoints (tread -> talus, talus -> cliff) for talus F-slope r.
    #[inline]
    pub fn breaks(&self, r: f64) -> (f64, f64) { let ut = self.ft / self.p; (ut, ut + if self.fb > 0.0 { self.fb / r } else { 0.0 }) }
}

/// The re-profile function of one bare-rock T1b post.
#[derive(Clone, Debug)]
pub struct Corner {
    pub r: f64,
    pub q: f64,
    /// Ascending period boundaries (m).
    pub ladder: Vec<f64>,
    /// Lane SHARP2: Some = the 3-limb (tread, talus, cliff) profile; None = ct_sharp_v1's 2-limb profile.
    pub t3: Option<Tread3>,
    /// Lane PERIOD-ENDS-2: Some = the published exact F (knot fractions, apron / massive period rules); None = r, q, t3 as before.
    pub ex: Option<Exact>,
}

/// Lane PERIOD-ENDS-2: a T1b post's exact F from pf-region ct_ends: knot fractions u = (h - s0) / T, f = (F - s0) / T at the
/// tread -> talus and talus -> cliff breaks (u[0] = f[0] = 0: no tread), the same in every period of the ladder.
#[derive(Clone, Copy, Debug)]
pub struct Exact {
    pub u: [f64; 2],
    pub f: [f64; 2],
    /// Apron (kind 5): the own period [s0, s1] (s0 = the hillslope base) overrides the ladder's [s_dn, s1] above s0; NaN otherwise.
    pub apron: (f64, f64),
    /// Massive (kind 6): F is the identity outside the ladder's [s0, s1].
    pub massive: bool,
}

impl Corner {
    /// Lane PERIOD-ENDS-2: a massive corner (kind 6) at a height outside its [s0, s1] (F is the identity there).
    #[inline]
    pub fn massive_out(&self, x: f64) -> bool { self.ex.map(|e| e.massive && (x <= self.ladder[0] || x >= self.ladder[self.ladder.len() - 1])).unwrap_or(false) }
    /// Lane PERIOD-ENDS-2: the corner has a tread limb.
    #[inline]
    pub fn has_tread(&self) -> bool { self.t3.is_some() || self.ex.map(|e| e.u[0] > 0.0).unwrap_or(false) }
    /// Lane PERIOD-ENDS-2: the F-slopes (tread or NaN, talus, cliff) of the corner's limbs from its knots.
    pub fn limb_slopes(&self) -> (f64, f64, f64) {
        let (u, f) = self.knots();
        let sl = |du: f64, df: f64| if du > 1e-9 { df / du } else { f64::NAN };
        (if u[0] > 1e-9 { f[0] / u[0] } else { f64::NAN }, sl(u[1] - u[0], f[1] - f[0]), sl(1.0 - u[1], 1.0 - f[1]))
    }
    /// Lane PERIOD-ENDS-2: the F-slope dF/dh at height x (1 outside a massive period).
    pub fn f_slope(&self, step: f64, x: f64) -> f64 {
        if self.massive_out(x) { return 1.0; }
        let (s0, s1) = self.period(step, x);
        let t = s1 - s0;
        if !(t > 0.0) { return 1.0; }
        let u = ((x - s0) / t).clamp(0.0, 1.0);
        let (uk, fk) = self.knots();
        let sl = |a: f64, b: f64, fa: f64, fb: f64| if b - a > 1e-9 { (fb - fa) / (b - a) } else { 1.0 };
        if u < uk[0] { sl(0.0, uk[0], 0.0, fk[0]) } else if u < uk[1] { sl(uk[0], uk[1], fk[0], fk[1]) } else { sl(uk[1], 1.0, fk[1], 1.0) }
    }
    #[inline]
    fn split(&self) -> (f64, f64) {
        let cc = (1.0 - 1.0 / self.r) / (1.0 / self.q - 1.0 / self.r);
        (cc, (1.0 - cc) / self.r)
    }
    /// The period [s0, s1] of the ladder (extended by the regular step period beyond its ends) holding height x.
    #[inline]
    pub fn period(&self, step: f64, x: f64) -> (f64, f64) {
        if let Some(e) = &self.ex { if x >= e.apron.0 && x < e.apron.1 { return e.apron; } }
        let l = &self.ladder;
        let k = l.partition_point(|v| *v <= x);
        if k == 0 {
            let s1 = l[0] - step * ((l[0] - x) / step).floor();
            (s1 - step, s1)
        } else if k == l.len() {
            let s0 = l[k - 1] + step * ((x - l[k - 1]) / step).floor();
            (s0, s0 + step)
        } else { (l[k - 1], l[k]) }
    }
    /// Lane PHASE: the knot fractions of F over one period: u = (h - s0) / T at the tread -> talus and talus -> cliff breaks, and
    /// (z - s0) / T there (2-limb: no tread, (0, ut) and (0, 1 - cc)).
    #[inline]
    pub fn knots(&self) -> ([f64; 2], [f64; 2]) {
        if let Some(e) = &self.ex { return (e.u, e.f); }
        if let Some(t3) = &self.t3 { let (ut, ub) = t3.breaks(self.r); ([ut, ub], [t3.ft, t3.ft + t3.fb]) } else { let (cc, ut) = self.split(); ([0.0, ut], [0.0, 1.0 - cc]) }
    }
    /// Lane TREADSOIL: the tread weight e at height x: 1 on the tread limb, falling linearly to 0 across the talus limb, 0 on
    /// the cliff limb; 0 for 2-limb corners.
    #[inline]
    pub fn tread_e(&self, step: f64, x: f64) -> f64 {
        if let Some(e) = &self.ex {
            if !(e.u[0] > 0.0) { return 0.0; }
            let (s0, s1) = self.period(step, x);
            let t = s1 - s0;
            if !(t > 0.0) || self.massive_out(x) { return 0.0; }
            return edge_e(((x - s0) / t).clamp(0.0, 1.0), e.u[0], e.u[1]);
        }
        let Some(t3) = &self.t3 else { return 0.0 };
        let (s0, s1) = self.period(step, x);
        let t = s1 - s0;
        if !(t > 0.0) { return 0.0; }
        let (ut, ub) = t3.breaks(self.r);
        edge_e(((x - s0) / t).clamp(0.0, 1.0), ut, ub)
    }
    /// Lane CELL-RISER: the inverse of F (exact corners, not apron / massive): the height h with F(h) = z on the ladder (F maps
    /// each period onto itself, identity at the levels, strictly increasing).
    pub fn f_inv(&self, step: f64, z: f64) -> f64 {
        let (s0, s1) = self.period(step, z);
        let t = s1 - s0;
        if !(t > 0.0) { return z; }
        let v = ((z - s0) / t).clamp(0.0, 1.0);
        let (uk, fk) = self.knots();
        let u = if v < fk[0] && fk[0] > 1e-12 { v * uk[0] / fk[0] }
            else if v < fk[1] && fk[1] - fk[0] > 1e-12 { uk[0] + (v - fk[0]) * (uk[1] - uk[0]) / (fk[1] - fk[0]) }
            else if 1.0 - fk[1] > 1e-12 { uk[1] + (v - fk[1]) * (1.0 - uk[1]) / (1.0 - fk[1]) } else { uk[1] };
        s0 + t * u
    }
    /// T1b's remap at height x on this ladder: (F(x), regime 1 talus / 2 cliff).
    pub fn f(&self, step: f64, x: f64) -> (f64, u8) {
        let (s0, s1) = self.period(step, x);
        if let Some(e) = &self.ex {
            if self.massive_out(x) { return (x, if x <= self.ladder[0] { 1 } else { 2 }); }
            let t = s1 - s0;
            if !(t > 0.0) { return (x, 1); }
            let u = ((x - s0) / t).clamp(0.0, 1.0);
            let (uk, fk) = (e.u, e.f);
            return if u < uk[0] { (s0 + t * u * fk[0] / uk[0], 3) }
                else if u < uk[1] { (s0 + t * (fk[0] + (u - uk[0]) * (fk[1] - fk[0]) / (uk[1] - uk[0])), 1) }
                else if 1.0 - uk[1] > 1e-12 { (s0 + t * (fk[1] + (u - uk[1]) * (1.0 - fk[1]) / (1.0 - uk[1])).min(1.0), 2) } else { (s1, 2) };
        }
        let t = s1 - s0;
        if !(t > 0.0) { return (x, 1); }
        if let Some(t3) = &self.t3 {
            let (ut, ub) = t3.breaks(self.r);
            let u = ((x - s0) / t).clamp(0.0, 1.0);
            return if u < ut { (s0 + t * u * t3.p, 3) } else if u < ub { (s0 + t * (t3.ft + (u - ut) * self.r), 1) } else { (s0 + t * ((t3.ft + t3.fb) + (u - ub) * self.q).min(1.0), 2) };
        }
        let (cc, ut) = self.split();
        let u = ((x - s0) / t).clamp(0.0, 1.0);
        if u < ut { (s0 + t * u * self.r, 1) } else { (s0 + t * ((1.0 - cc) + (u - ut) * self.q).min(1.0), 2) }
    }
}

/// The 2 m re-profile over a T1b window lattice.
pub struct SharpField {
    /// Lane SPLIT-2M (ct_sharp_blend_period=1): knots blend only between corners on the same period (SHARP_BLEND_PERIOD_RULE).
    pub same_period: bool,
    pub n: usize,
    pub gx0: i64,
    pub gy0: i64,
    pub step_m: f64,
    /// T1b staircase shift z - h per post (m).
    pub d: Vec<f64>,
    /// Index into `corners` per post (u32::MAX: identity).
    pub corner: Vec<u32>,
    pub corners: Vec<Corner>,
    /// Lane TREAD: the cliff tangent ta (rock-mass multiplier applied) per corner (parallel to `corners`).
    pub corner_ta: Vec<f64>,
    /// Lane SHARP2: the window carries pf-region ct_tread's fields (3-limb corners possible; 2 m regime 3 = tread).
    pub tread: bool,
    pub stats: serde_json::Value,
    /// Lane PHASE: 0 = F-value blend (default), 1 = ct_sharp_blend=phase, 2 = ct_sharp_blend=nearest (diagnostic).
    pub blend: u8,
    /// Lane PHASE: the merge tolerance (m): corners whose period boundaries agree within it share one blended F.
    pub tol: f64,
    /// Lane TREADSOIL: the mantle edge e0 (soil weight where the colluvium wedge starts on a tread limb; 0 = off, the linear fade).
    pub soil_edge: f64,
    /// Lane TREADSOIL: the wedge's extra-gradient limit (0 = none).
    pub soil_tw: f64,
    /// Lane TREADSOIL: T1b regime 0 (soil-mantled) per post (empty when soil_edge is off).
    pub soil: Vec<bool>,
    /// Lane PERIOD-ENDS-2 (ct_sharp_ends_label=1): the 2 m regime is the limb of the blended F at h2 (SHARP_ENDS_LABEL_RULE).
    pub ends_label: bool,
    /// Lane CELL-RISER: the pseudo-base shift dh per T1b post (n x n, 0 where none; empty when ct_sharp_cell_riser is off).
    pub dh_post: Vec<f64>,
    /// Lane CELL-RISER: dh interpolation, 0 bilinear, 1 Catmull-Rom bicubic, 2 smoothstep-weighted bilinear (C1, no overshoot), 3 monotone cubic (CELL-RISER-2).
    pub dh_mode: u8,
}

/// Lane PERIOD-ENDS-2 (SHARP_ENDS_LABEL_RULE): the limb of a blended F-slope `sb` against the reference corner's limb F-slopes
/// (tread p or NaN, talus r, cliff q): cliff above the geometric mean of r and q, tread below that of p and r, else talus;
/// `fallback` where the reference has no talus or cliff limb.
#[inline]
pub fn limb_of_slope(sb: f64, (p, r, q): (f64, f64, f64), fallback: u8) -> u8 {
    if !(r.is_finite() && q.is_finite() && r > 0.0 && q > 0.0) { return fallback; }
    if sb >= (r * q).sqrt() { 2 } else if p.is_finite() && p > 0.0 && sb <= (p * r).sqrt() { 3 } else { 1 }
}

#[inline]
fn cr(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [(-t3 + 2.0 * t2 - t) * 0.5, (3.0 * t3 - 5.0 * t2 + 2.0) * 0.5, (-3.0 * t3 + 4.0 * t2 + t) * 0.5, (t3 - t2) * 0.5]
}

/// Lane TREADSOIL: the mantle edge (soil weight) where the colluvium wedge starts on a tread limb.
pub const SOIL_EDGE_E0: f64 = 0.5;

/// Lane CELL-RISER: default D (m): a post's cliff foot more than this above the lowest same-period neighbour foot is re-profiled.
pub const CELL_RISER_D_M: f64 = 20.0;
/// Lane CELL-RISER: the default dh interpolation (ct_sharp_cell_riser_dh).
pub const CELL_RISER_DH_DEFAULT: &str = "smooth";
/// Lane CELL-RISER: dh interpolation code of a ct_sharp_cell_riser_dh value (None = the default).
pub fn cell_riser_dh_mode(v: Option<&str>) -> u8 { match v.unwrap_or(CELL_RISER_DH_DEFAULT) { "bilinear" => 0, "bicubic" => 1, "monotone" => 3, _ => 2 } }

/// Lane CELL-RISER-2: monotone cubic Hermite through v[1] (t = 0) and v[2] (t = 1) on unit spacing, Fritsch-Butland
/// harmonic-mean tangents (0 at a local extremum or a flat interval: no overshoot, monotone where the data are).
#[inline]
pub fn mono_cubic(v: [f64; 4], t: f64) -> f64 {
    let hm = |a: f64, b: f64| if a * b <= 0.0 { 0.0 } else { 2.0 * a * b / (a + b) };
    let d0 = v[2] - v[1];
    let (m1, m2) = (hm(v[1] - v[0], d0), hm(d0, v[3] - v[2]));
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * v[1] + (t3 - 2.0 * t2 + t) * m1 + (-2.0 * t3 + 3.0 * t2) * v[2] + (t3 - t2) * m2
}

/// Lane STEEP-TALUS: default threshold angle (deg) of ct_sharp_talus_rock: steeper than this over 8 m both up- and downslope a
/// regime-1 post is rock (coarse talus rests at 34-40 deg; the DE talus repose is 0.75 = 36.9 deg; 45 deg leaves the 2 m
/// staircase's own talus spread, r |grad h2| p90 ~ 44 deg, to the talus).
pub const TALUS_ROCK_DEG: f64 = 45.0;
/// Lane STEEP-TALUS: the one-sided baseline of the ct_sharp_talus_rock slope (posts; 4 posts = 8 m).
pub const TALUS_ROCK_POSTS: usize = 4;

/// Lane STEEP-TALUS (SHARP_TALUS_ROCK_RULE): relabel the regime-1 posts of an n x n domain surface `z` (2 m posts) whose minmod
/// 8 m slope exceeds `tan_thr` as regime 2. Per axis the minmod of the one-sided TALUS_ROCK_POSTS differences (0 when they
/// disagree in sign); posts closer than TALUS_ROCK_POSTS to the ring are not tested. Returns the relabelled posts.
pub fn talus_rock_relabel(n: usize, z: &[f64], rg: &mut [u8], tan_thr: f64) -> Vec<usize> {
    let m = TALUS_ROCK_POSTS;
    let base = m as f64 * crate::POST_M;
    let mm = |a: f64, b: f64| if a * b > 0.0 { a.signum() * a.abs().min(b.abs()) } else { 0.0 };
    let mut out = Vec::new();
    if n <= 2 * m { return out; }
    for j in m..n - m { for i in m..n - m {
        let k = j * n + i;
        if rg[k] != 1 { continue; }
        let gx = mm((z[k] - z[k - m]) / base, (z[k + m] - z[k]) / base);
        let gy = mm((z[k] - z[k - m * n]) / base, (z[k + m * n] - z[k]) / base);
        if gx.hypot(gy) > tan_thr { out.push(k); }
    } }
    for &k in &out { rg[k] = 2; }
    out
}
/// Lane TREADSOIL: the largest extra gradient (tangent) of the colluvium wedge between the mantle edge and the soil corner;
/// a thicker wedge moves the mantle edge back toward the soil corner's linear fade.
pub const SOIL_EDGE_TW: f64 = 1.0;

/// Lane TREADSOIL: tread weight at position x between the tread -> talus break a and the talus -> cliff break b.
#[inline]
fn edge_e(x: f64, a: f64, b: f64) -> f64 { if x < a { 1.0 } else if x < b && b - a > 1e-12 { (b - x) / (b - a) } else { 0.0 } }

impl SharpField {
    /// Lane TREADSOIL: the factor on corner k's w_k (Fk - h2) at soil weight ws with tread weight e (1 when off or ws = 0).
    /// `dz` = Fk - h2 (the cut of corner k at the post): the mantle edge e0 shrinks with the wedge it would carry,
    /// e0 = min(SOIL_EDGE_E0, 1 - |dz| / (T1B_M TW)), so the wedge's extra gradient |dz| / ((1 - e0) 32 m) stays within
    /// TW (soil_tw, default SOIL_EDGE_TW; 0 = no limit) (e0 = 0 at |dz| >= 32 m x TW: the linear fade).
    #[inline]
    fn edge_g(&self, ws: f64, e: f64, dz: f64) -> f64 {
        if self.soil_edge <= 0.0 || ws <= 0.0 || e <= 0.0 { return 1.0; }
        let e0 = if self.soil_tw > 0.0 { self.soil_edge.min(1.0 - dz.abs() / (crate::T1B_M * self.soil_tw)) } else { self.soil_edge };
        if e0 <= 0.0 { return 1.0; }
        (1.0 / (1.0 - ws).max(1e-9)).min(1.0 / (1.0 - e0 * e))
    }

    pub fn build(c: &CtInputs, p: &SharpParams) -> SharpField {
        let n = c.n;
        let nn = n * n;
        let dx = crate::T1B_M;
        let at = |x: i64, y: i64| c.h[y as usize * n + x as usize];
        let d: Vec<f64> = (0..nn).map(|i| c.z[i] - c.h[i]).collect();
        // per post: (r, q, own boundary, regime) for bare rock
        let mut tent: Vec<Option<(f64, f64, f64, u8)>> = vec![None; nn];
        let mut tant = vec![0.0f64; nn];
        let (mut n_bare, mut n_fixg, mut n_steep, mut n_flat) = (0u64, 0u64, 0u64, 0u64);
        // lane SHARP2: per post the 3-limb profile (None: 2-limb, every step as ct_sharp_v1)
        let mut t3v: Vec<Option<Tread3>> = vec![None; nn];
        let (mut n3, mut n3_fixg, mut n3_flat, mut n3_rg3_noft, mut n3_talus_dpos) = (0u64, 0u64, 0u64, 0u64, 0u64);
        let mut n3_rg = [0u64; 4];
        for i in 0..nn {
            let rg = c.regime[i];
            if rg == 0 { continue; }
            let cf = c.cliff_permille[i] as f64 / 1000.0;
            if cf >= 1.0 { n_steep += 1; continue; }
            if let Some((ftp, ttp)) = c.tread {
                if ftp[i] > 0 || rg == 3 {
                    // lane SHARP2: 3-limb post
                    let code = c.litho[i] as usize;
                    let Some(k) = c.classes.get(code).copied() else { continue };
                    if k.mode == 0 || !(k.tb > 0.0) { continue; }
                    let (ft0, tt) = (ftp[i] as f64 / 1000.0, ttp[i] as f64 / 1000.0);
                    if !(ft0 > 0.0 && tt > 0.0 && tt < k.tb) { n3_rg3_noft += 1; continue; }
                    let fb = (1.0 - ft0 - cf).max(0.0);
                    let ft = 1.0 - cf - fb;
                    let (x, y) = ((i % n) as i64, (i / n) as i64);
                    let d1 = |a: Option<f64>, b: Option<f64>, m: f64| match (a, b) { (Some(p), Some(q)) => (q - p) / (2.0 * dx), (None, Some(q)) => (q - m) / dx, (Some(p), None) => (m - p) / dx, _ => 0.0 };
                    let m = c.h[i];
                    let gx = d1(if x > 0 { Some(at(x - 1, y)) } else { None }, if x + 1 < n as i64 { Some(at(x + 1, y)) } else { None }, m);
                    let gy = d1(if y > 0 { Some(at(x, y - 1)) } else { None }, if y + 1 < n as i64 { Some(at(x, y + 1)) } else { None }, m);
                    let mut g = gx.hypot(gy);
                    let tb = k.tb;
                    let mut ta = cf / (1.0 / g - ft / tt - fb / tb);
                    let rm = ta / k.ta;
                    if !(ta.is_finite() && ta > tb && rm >= c.rockmass.0 - 1e-9 && rm <= c.rockmass.1 + 1e-9) {
                        ta = (k.ta * if rm.is_finite() && rm > 0.0 { rm.clamp(c.rockmass.0, c.rockmass.1) } else { 1.0 }).max(tb * 1.0001);
                        g = 1.0 / (ft / tt + fb / tb + cf / ta);
                        n3_fixg += 1;
                    }
                    let (p, r, q) = (tt / g, tb / g, ta / g);
                    if !(p < r && p < 1.0 && q > 1.0) { n3_flat += 1; continue; }
                    let t3 = Tread3 { p, ft, fb, fc: cf };
                    // own boundary: tread -> s0 (exact), cliff -> s1 (exact), talus -> NaN (solved from the ladder below)
                    let own = match rg { 3 => m + d[i] / (1.0 - p), 2 => m + d[i] / (1.0 - q), _ => f64::NAN };
                    if rg == 1 && !(d[i] < 0.0) { n3_talus_dpos += 1; continue; }
                    tent[i] = Some((r, q, own, rg));
                    tant[i] = ta;
                    t3v[i] = Some(t3);
                    n3 += 1;
                    n3_rg[rg.min(3) as usize] += 1;
                    n_bare += 1;
                    continue;
                }
            }
            let code = c.litho[i] as usize;
            let Some(k) = c.classes.get(code).copied() else { continue };
            if k.mode == 0 || !(k.tb > 0.0) { continue; }
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            let d1 = |a: Option<f64>, b: Option<f64>, m: f64| match (a, b) { (Some(p), Some(q)) => (q - p) / (2.0 * dx), (None, Some(q)) => (q - m) / dx, (Some(p), None) => (m - p) / dx, _ => 0.0 };
            let m = c.h[i];
            let gx = d1(if x > 0 { Some(at(x - 1, y)) } else { None }, if x + 1 < n as i64 { Some(at(x + 1, y)) } else { None }, m);
            let gy = d1(if y > 0 { Some(at(x, y - 1)) } else { None }, if y + 1 < n as i64 { Some(at(x, y + 1)) } else { None }, m);
            let mut g = gx.hypot(gy);
            let tb = k.tb;
            let mut ta = cf / (1.0 / g - (1.0 - cf) / tb);
            let rm = ta / k.ta;
            if !(ta.is_finite() && ta > tb && rm >= c.rockmass.0 - 1e-9 && rm <= c.rockmass.1 + 1e-9) {
                ta = (k.ta * if rm.is_finite() { rm.clamp(c.rockmass.0, c.rockmass.1) } else { 1.0 }).max(tb * 1.0001);
                g = 1.0 / (cf / ta + (1.0 - cf) / tb);
                n_fixg += 1;
            }
            let (r, q) = (tb / g, ta / g);
            if !(r < 1.0 && q > 1.0) { n_flat += 1; continue; }
            let mm = if rg == 1 { r } else { q };
            tent[i] = Some((r, q, m + d[i] / (1.0 - mm), rg));
            tant[i] = ta;
            n_bare += 1;
        }
        let rad = p.radius_posts;
        let tol = p.merge_tol_m;
        let step = c.step_m;
        let mut corner = vec![u32::MAX; nn];
        let mut corners: Vec<Corner> = Vec::new();
        let mut corner_ta: Vec<f64> = Vec::new();
        let (mut n_ext, mut n_drop, mut lev_sum) = (0u64, 0u64, 0u64);
        let mut cand: Vec<(f64, bool)> = Vec::new();
        // lane SHARP2 ladder counters: talus period from a level above / below / the step period / centred, failed, dropped, extended
        let mut t3stat = [0u64; 7];
        for i in 0..nn {
            let Some((r, q, own, rg)) = tent[i] else { continue };
            // the post's gradient G = tb / r (the class talus tangent over the talus-limb slope of F)
            let tbk = c.classes.get(c.litho[i] as usize).map(|k| k.tb).unwrap_or(0.7);
            let tb_of = |r: f64| tbk / r;
            let (x, y) = ((i % n) as i64, (i / n) as i64);
            cand.clear();
            for yy in (y - rad).max(0)..=(y + rad).min(n as i64 - 1) {
                for xx in (x - rad).max(0)..=(x + rad).min(n as i64 - 1) {
                    let j = yy as usize * n + xx as usize;
                    if let Some(t) = tent[j] { if t.2.is_finite() { cand.push((t.2, j == i)); } }
                }
            }
            cand.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            if let Some(t3) = t3v[i] {
                // lane SHARP2: the 3-limb ladder
                let Some(lad) = ladder3(&cand, c.h[i], d[i], own, rg, r, t3, step, tol, rad as f64 * dx * tb_of(r), p.t3_consensus, &mut t3stat) else { t3stat[4] += 1; continue };
                lev_sum += lad.len() as u64;
                corner[i] = corners.len() as u32;
                corners.push(Corner { r, q, ladder: lad, t3: Some(t3), ex: None });
                corner_ta.push(tant[i]);
                continue;
            }
            // merge into levels; the own boundary's cluster is replaced by the own value (exact)
            let mut lad: Vec<f64> = Vec::new();
            let mut s = 0usize;
            while s < cand.len() {
                let mut e = s + 1;
                while e < cand.len() && cand[e].0 - cand[e - 1].0 < tol { e += 1; }
                if !cand[s..e].iter().any(|v| v.1) {
                    let m = e - s;
                    let med = if m % 2 == 1 { cand[s + m / 2].0 } else { 0.5 * (cand[s + m / 2 - 1].0 + cand[s + m / 2].0) };
                    lad.push(med);
                }
                s = e;
            }
            let h = c.h[i];
            let (lo, hi) = (own.min(h), own.max(h));
            lad.retain(|v| !(lo - tol < *v && *v < hi + tol));
            lad.push(own);
            lad.sort_by(|a, b| a.total_cmp(b));
            let cn = Corner { r, q, ladder: Vec::new(), t3: None, ex: None };
            let (cc, ut) = cn.split();
            if rg == 1 {
                loop {
                    let Some(&up) = lad.iter().find(|v| **v > own) else { break };
                    if h <= own + ut * (up - own) + 1e-9 { break; }
                    lad.retain(|v| *v != up); n_drop += 1;
                }
                // no boundary above within the radius: the riser foot is at least `radius` posts upslope (no bare neighbour there
                // shows a riser), at least one regular period
                let reach = rad as f64 * dx * tb_of(r);
                if !lad.iter().any(|v| *v > own) { lad.push(own + step.max((h - own + reach) * r / (1.0 - cc))); n_ext += 1; }
            } else {
                loop {
                    let Some(&dn) = lad.iter().rev().find(|v| **v < own) else { break };
                    if h >= dn + ut * (own - dn) - 1e-9 { break; }
                    lad.retain(|v| *v != dn); n_drop += 1;
                }
                let reach = rad as f64 * dx * tb_of(r);
                if !lad.iter().any(|v| *v < own) { lad.insert(0, own - step.max((own - h + reach) / (1.0 - ut))); n_ext += 1; }
            }
            lev_sum += lad.len() as u64;
            corner[i] = corners.len() as u32;
            corners.push(Corner { r, q, ladder: lad, t3: None, ex: None });
            corner_ta.push(tant[i]);
        }
        let stats = serde_json::json!({"bare_posts": n_bare, "ta_clamped_posts": n_fixg, "steeper_than_cliff_posts": n_steep, "not_steeper_than_talus_posts": n_flat,
            "ladder_levels_mean": if n_bare > 0 { lev_sum as f64 / n_bare as f64 } else { 0.0 }, "levels_dropped_inconsistent": n_drop, "regular_period_extensions": n_ext});
        let mut stats = stats;
        if c.tread.is_some() {
            stats["tread"] = serde_json::json!({"rule": SHARP_TREAD_RULE, "posts_3limb": n3, "posts_3limb_by_t1b_regime": n3_rg, "ta_clamped_posts": n3_fixg,
                "not_ordered_posts": n3_flat, "tread_field_missing_posts": n3_rg3_noft, "talus_not_below_h_posts": n3_talus_dpos,
                "talus_period_from_level_above": t3stat[0], "talus_period_from_level_below": t3stat[1], "talus_period_step": t3stat[2], "talus_period_centred": t3stat[3],
                "ladder_failed_posts": t3stat[4], "levels_dropped_inconsistent": t3stat[5], "regular_period_extensions": t3stat[6]});
        }
        // lane PERIOD-ENDS-2 (ct_sharp_ends=1): every bare T1b post with published period ends and knots (not flagged held /
        // moved, kind + 128) takes its exact F; the other posts keep the corner built above (the neighbour candidates `tent` are
        // unchanged, so their ladders are as without the key)
        if let (Some(_), Some(e)) = (p.ends.as_ref(), c.ends) {
            let none = i32::MIN;
            let mm = |v: i32| v as f64 / 1000.0;
            let mut st = [0u64; 10];
            for i in 0..nn {
                if c.regime[i] == 0 { continue; }
                let kind = e.kind[i];
                if e.s0[i] == none || e.s1[i] == none || kind & 0x7f == 0 { st[1] += 1; continue; }
                if kind & 0x80 != 0 { st[2] += 1; continue; }
                if [e.h1[i], e.z1[i], e.h2[i], e.z2[i]].contains(&none) { st[1] += 1; continue; }
                let (s0, s1) = (mm(e.s0[i]), mm(e.s1[i]));
                let t = s1 - s0;
                if !(t > 0.0) { st[1] += 1; continue; }
                let u = [((mm(e.h1[i]) - s0) / t).clamp(0.0, 1.0), ((mm(e.h2[i]) - s0) / t).clamp(0.0, 1.0)];
                let f = [((mm(e.z1[i]) - s0) / t).clamp(0.0, 1.0), ((mm(e.z2[i]) - s0) / t).clamp(0.0, 1.0)];
                let (u, f) = if u[0] * t < 1e-3 { ([0.0, u[1]], [0.0, f[1]]) } else { (u, f) };
                let (dn, up) = (e.s_dn[i], e.s_up[i]);
                if dn == none { st[3] += 1; }
                if up == none { st[4] += 1; }
                let k = kind & 0x7f;
                let massive = k == 6;
                let mut lad: Vec<f64> = Vec::with_capacity(4);
                if dn != none && !massive && mm(dn) < s0 { lad.push(mm(dn)); }
                let apron = if k == 5 { st[5] += 1; (s0, s1) } else { lad.push(s0); (f64::NAN, f64::NAN) };
                lad.push(s1);
                if up != none && !massive && mm(up) > s1 { lad.push(mm(up)); }
                lad.sort_by(|a, b| a.total_cmp(b));
                lad.dedup();
                if massive { st[6] += 1; }
                let ex = Exact { u, f, apron, massive };
                // cliff tangent: ta = S x the cliff limb's F-slope (1 - f2) / (1 - u2); else the old corner's / the class tangent
                let sg = e.grad_permille[i] as f64 / 1000.0;
                let q = if 1.0 - u[1] > 1e-6 { (1.0 - f[1]) / (1.0 - u[1]) } else { f64::NAN };
                let ta_k = c.classes.get(c.litho[i] as usize).map(|k| k.ta).unwrap_or(0.0);
                let ta = if sg > 0.0 && q.is_finite() && q > 1.0 { sg * q } else { st[7] += 1; if corner[i] != u32::MAX { corner_ta[corner[i] as usize] } else { ta_k } };
                let r = if u[1] - u[0] > 1e-9 { (f[1] - f[0]) / (u[1] - u[0]) } else { 1.0 };
                let cn = Corner { r, q: if q.is_finite() { q } else { 1.0 }, ladder: lad, t3: None, ex: Some(ex) };
                if corner[i] == u32::MAX { st[8] += 1; corner[i] = corners.len() as u32; corners.push(cn); corner_ta.push(ta); }
                else { let ci = corner[i] as usize; corners[ci] = cn; corner_ta[ci] = ta; }
                st[0] += 1;
                if u[0] > 0.0 { st[9] += 1; }
            }
            stats["ends"] = serde_json::json!({"rule": SHARP_ENDS_RULE, "exact_posts": st[0], "exact_posts_3limb": st[9], "exact_posts_new_corner": st[8],
                "bare_posts_without_ends": st[1], "held_posts_old_path": st[2], "s_dn_missing_posts": st[3], "s_up_missing_posts": st[4],
                "apron_posts": st[5], "massive_posts": st[6], "ta_fallback_posts": st[7]});
        }
        // lane CELL-RISER (ct_sharp_cell_riser=1): same-period posts whose cliff foot sits more than D above the lowest foot of
        // their 8-neighbours take that neighbour's knots and keep their height by a pseudo-base shift dh (SHARP_CELL_RISER_RULE)
        let mut dh_post: Vec<f64> = Vec::new();
        if let (Some(_), Some(e), Some(_)) = (p.cell_riser.as_ref(), c.ends, p.ends.as_ref()) {
            let dthr = p.cell_riser_m.unwrap_or(CELL_RISER_D_M);
            dh_post = vec![0.0; nn];
            let mm = |v: i32| v as f64 / 1000.0;
            // eligible: an exact corner on a plain period (not apron, not massive): (s0, s1, cliff foot z2)
            let el: Vec<Option<(f64, f64, f64)>> = (0..nn).map(|i| {
                let ci = corner[i];
                if ci == u32::MAX { return None; }
                let ck = &corners[ci as usize];
                let ex = ck.ex?;
                if ex.massive || ex.apron.0.is_finite() { return None; }
                let (s0, s1) = (mm(e.s0[i]), mm(e.s1[i]));
                Some((s0, s1, s0 + (s1 - s0) * ex.f[1]))
            }).collect();
            let mut upd: Vec<(usize, usize)> = Vec::new();
            let (mut n_el, mut n_low) = (0u64, 0u64);
            let donor_mode: u8 = match p.cell_riser_donor.as_deref() { Some("median") => 1, Some("pair") => 2, _ => 0 };
            let mut n_unconf = 0u64;
            for i in 0..nn {
                let Some((s0, s1, z2)) = el[i] else { continue };
                n_el += 1;
                let (x, y) = ((i % n) as i64, (i / n) as i64);
                let mut best: Option<(f64, usize)> = None;
                let mut cand: Vec<(f64, usize)> = Vec::new();
                for yy in (y - 1).max(0)..=(y + 1).min(n as i64 - 1) { for xx in (x - 1).max(0)..=(x + 1).min(n as i64 - 1) {
                    let j = yy as usize * n + xx as usize;
                    if j == i { continue; }
                    let Some((t0, t1, zj)) = el[j] else { continue };
                    if (t0 - s0).abs() >= tol || (t1 - s1).abs() >= tol { continue; }
                    if donor_mode != 0 { cand.push((zj, j)); }
                    if best.map(|b| zj < b.0).unwrap_or(true) { best = Some((zj, j)); }
                } }
                // lane CELL-RISER-2: the donor by the median / the confirmed lowest foot (CELL_RISER2_RULE); stable order (z, index)
                if donor_mode != 0 {
                    cand.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                    best = if donor_mode == 1 { cand.get((cand.len().max(1) - 1) / 2).copied() }
                        else { (0..cand.len().saturating_sub(1)).find(|&k| cand[k + 1].0 - cand[k].0 <= dthr).map(|k| cand[k]) };
                    if best.is_none() && !cand.is_empty() { n_unconf += 1; }
                }
                let Some((zb, j)) = best else { continue };
                if z2 - zb > dthr { upd.push((i, j)); } else if zb >= z2 { n_low += 1; }
            }
            let mut dhs: Vec<f64> = Vec::with_capacity(upd.len());
            let mut n_cliff = 0u64;
            let mut n_capped = 0u64;
            let newc: Vec<(usize, usize, Corner, f64, f64)> = upd.iter().filter_map(|&(i, j)| {
                let (ci, cj) = (corner[i] as usize, corner[j] as usize);
                let (own, don) = (&corners[ci], &corners[cj]);
                let dex = don.ex.expect("eligible donor");
                let ex = Exact { u: dex.u, f: dex.f, apron: (f64::NAN, f64::NAN), massive: false };
                let cn = Corner { r: don.r, q: don.q, ladder: own.ladder.clone(), t3: None, ex: Some(ex) };
                let hh = cn.f_inv(step, c.z[i]);
                // lane CELL-RISER-2: the break moves |dh| / S horizontally; beyond cap x 32 m the post keeps its own knots
                if let Some(cap) = p.cell_riser_cap {
                    let sg = e.grad_permille[i] as f64 / 1000.0;
                    if !((hh - c.h[i]).abs() <= cap * crate::T1B_M * sg) { n_capped += 1; return None; }
                }
                if cn.f(step, hh).1 == 2 { n_cliff += 1; }
                Some((i, ci, cn, corner_ta[cj], hh - c.h[i]))
            }).collect();
            let n_re = newc.len();
            for (i, ci, cn, ta, dh) in newc { corners[ci] = cn; corner_ta[ci] = ta; dh_post[i] = dh; dhs.push(dh.abs()); }
            dhs.sort_by(|a, b| a.total_cmp(b));
            let qv = |f: f64| if dhs.is_empty() { 0.0 } else { dhs[((dhs.len() - 1) as f64 * f).round() as usize] };
            stats["cell_riser"] = serde_json::json!({"rule": SHARP_CELL_RISER_RULE, "d_m": dthr, "eligible_posts": n_el, "reprofiled_posts": n_re,
                "reprofiled_on_cliff_limb_posts": n_cliff, "lowest_foot_posts": n_low, "dh_abs_p50_m": qv(0.5), "dh_abs_p90_m": qv(0.9), "dh_abs_max_m": qv(1.0)});
            if p.cell_riser_donor.is_some() || p.cell_riser_cap.is_some() || p.cell_riser_dh.as_deref() == Some("monotone") {
                stats["cell_riser"]["v2"] = serde_json::json!({"rule": CELL_RISER2_RULE, "donor": p.cell_riser_donor.as_deref().unwrap_or("min"), "cap_cells": p.cell_riser_cap,
                    "candidates_over_d": upd.len(), "capped_posts": n_capped, "unconfirmed_posts": n_unconf});
            }
        }
        let blend = match p.blend.as_deref() { None => 0, Some(r) if r == SHARP_BLEND_PHASE_RULE => 1, Some(_) => 2 };
        let soil_edge = if p.soil_edge.is_some() { SOIL_EDGE_E0 } else { 0.0 };
        let soil = if soil_edge > 0.0 { c.regime.iter().map(|&r| r == 0).collect() } else { Vec::new() };
        SharpField { same_period: p.blend_period.is_some(), n, gx0: c.gx0, gy0: c.gy0, step_m: step, d, corner, corners, corner_ta, tread: c.tread.is_some(), stats, blend, tol, soil_edge, soil_tw: p.soil_edge_tw.unwrap_or(SOIL_EDGE_TW), soil, ends_label: p.ends_label.is_some(), dh_post, dh_mode: cell_riser_dh_mode(p.cell_riser_dh.as_deref()) }
    }

    /// Lane TREADSOIL: soil corners on or next to (8-neighbour) a fixed T1b post (channel corridor, sea, other fixed posts:
    /// valley-floor and shore soil, not bench colluvium) keep the linear fade, so an extended bench never undercuts a channel corridor.
    pub fn exclude_soil_near_fixed(&mut self, fixed: &[u8]) {
        if self.soil.is_empty() || fixed.len() != self.soil.len() { return; }
        let n = self.n as i64;
        let mut cnt = 0u64;
        for y in 0..n { for x in 0..n {
            let i = (y * n + x) as usize;
            if !self.soil[i] { continue; }
            let near = (-1..=1).any(|dy| (-1..=1).any(|dx| { let (xx, yy) = (x + dx, y + dy); xx >= 0 && yy >= 0 && xx < n && yy < n && fixed[(yy * n + xx) as usize] != 0 }));
            if near { self.soil[i] = false; cnt += 1; }
        } }
        self.stats["soil_edge_posts"] = serde_json::json!(self.soil.iter().filter(|v| **v).count());
        self.stats["soil_edge_excluded_near_fixed_posts"] = serde_json::json!(cnt);
    }

    #[inline]
    fn clampi(&self, gx: i64, gy: i64) -> usize {
        let x = (gx - self.gx0).clamp(0, self.n as i64 - 1) as usize;
        let y = (gy - self.gy0).clamp(0, self.n as i64 - 1) as usize;
        y * self.n + x
    }

    /// Bicubic (Catmull-Rom, as B1) of the T1b shift field at chart (x, y) m.
    pub fn shift_bicubic(&self, x: f64, y: f64) -> f64 { self.bicubic_of(&self.d, x, y) }

    /// Lane CELL-RISER: the pseudo-base shift at chart (x, y) m from the per-post dh (dh_mode: 0 bilinear, 1 Catmull-Rom bicubic,
    /// 2 bilinear with smoothstep weights 3t^2 - 2t^3: passes through every post, C1 across cell edges, no overshoot; 3 monotone cubic, lane CELL-RISER-2).
    pub fn dh_at(&self, x: f64, y: f64) -> f64 {
        if self.dh_mode == 1 { return self.bicubic_of(&self.dh_post, x, y); }
        let (fx, fy) = (x / crate::T1B_M, y / crate::T1B_M);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        if self.dh_mode == 3 {
            let (tx, ty) = (fx - ix as f64, fy - iy as f64);
            let row = |gy: i64| mono_cubic([0, 1, 2, 3].map(|q| self.dh_post[self.clampi(ix - 1 + q, gy)]), tx);
            return mono_cubic([row(iy - 1), row(iy), row(iy + 1), row(iy + 2)], ty);
        }
        let (mut tx, mut ty) = (fx - ix as f64, fy - iy as f64);
        if self.dh_mode == 2 { tx = tx * tx * (3.0 - 2.0 * tx); ty = ty * ty * (3.0 - 2.0 * ty); }
        let v = |gx: i64, gy: i64| self.dh_post[self.clampi(gx, gy)];
        (1.0 - ty) * ((1.0 - tx) * v(ix, iy) + tx * v(ix + 1, iy)) + ty * ((1.0 - tx) * v(ix, iy + 1) + tx * v(ix + 1, iy + 1))
    }

    /// Bicubic (Catmull-Rom, as B1) of a per-T1b-post field `v` (n x n) at chart (x, y) m (lane CELL-RISER: also the pseudo-base shift).
    pub fn bicubic_of(&self, v: &[f64], x: f64, y: f64) -> f64 {
        let (fx, fy) = (x / crate::T1B_M, y / crate::T1B_M);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (wx, wy) = (cr(fx - ix as f64), cr(fy - iy as f64));
        let mut acc = 0.0;
        for (j, wyj) in wy.iter().enumerate() {
            let mut row = 0.0;
            for (i, wxi) in wx.iter().enumerate() { row += wxi * v[self.clampi(ix - 1 + i as i64, iy - 1 + j as i64)]; }
            acc += wyj * row;
        }
        acc
    }

    /// The re-profiled height of a 2 m post at chart (x, y) m with pipeline height z2 (m): Some((new height, regime)), or None
    /// where the post is untouched (no bare-rock T1b post in the cell and no staircase shift in the bicubic stencil).
    pub fn apply(&self, x: f64, y: f64, z2: f64) -> Option<(f64, u8)> { self.apply_ta(x, y, z2).map(|v| (v.0, v.1)) }

    /// `apply` plus the cliff tangent at the post (bilinear over the cell's bare-rock T1b posts; 0 where none) (lane TREAD).
    pub fn apply_ta(&self, x: f64, y: f64, z2: f64) -> Option<(f64, u8, f64)> {
        let (fx, fy) = (x / crate::T1B_M, y / crate::T1B_M);
        let (ix, iy) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - ix as f64, fy - iy as f64);
        let cs = [(ix, iy, (1.0 - tx) * (1.0 - ty)), (ix + 1, iy, tx * (1.0 - ty)), (ix, iy + 1, (1.0 - tx) * ty), (ix + 1, iy + 1, tx * ty)];
        let s = self.shift_bicubic(x, y);
        let bare = cs.iter().any(|&(gx, gy, _)| self.corner[self.clampi(gx, gy)] != u32::MAX);
        if s == 0.0 && !bare { return None; }
        let h2 = z2 - s;
        // lane TREADSOIL: the bilinear weight of the cell's soil-mantled (T1b regime 0) corners (0 when the key is off)
        let ws: f64 = if self.soil_edge > 0.0 { cs.iter().filter(|c| c.2 > 0.0 && self.soil[self.clampi(c.0, c.1)]).map(|c| c.2).sum() } else { 0.0 };
        if self.blend != 0 {
            // lane PHASE: one F per post (phase) or the dominant corner's F (nearest)
            let mut bc: [(usize, f64); 4] = [(0, 0.0); 4];
            let (mut nb, mut wb, mut tsum, mut iref) = (0usize, 0.0, 0.0, 0usize);
            for &(gx, gy, wgt) in &cs {
                let k = self.corner[self.clampi(gx, gy)];
                if k == u32::MAX || wgt <= 0.0 { continue; }
                bc[nb] = (k as usize, wgt);
                if wgt > bc[iref].1 { iref = nb; }
                nb += 1;
                wb += wgt;
                tsum += wgt * self.corner_ta[k as usize];
            }
            let ta = if wb > 0.0 { tsum / wb } else { 0.0 };
            if nb == 0 { return Some((h2, 0, 0.0)); }
            // lane CELL-RISER: the corners' F are evaluated at the pseudo base h2 + dh (Catmull-Rom bicubic of the per-post shifts)
            let ht = if self.dh_post.is_empty() { h2 } else { h2 + self.dh_at(x, y) };
            if nb == 1 || self.blend == 2 {
                let ck = &self.corners[bc[iref].0];
                let (fz, rg) = ck.f(self.step_m, ht);
                let g = if ws > 0.0 && self.blend == 1 { self.edge_g(ws, ck.tread_e(self.step_m, ht), fz - h2) } else { 1.0 };
                return Some((h2 + g * wb * (fz - h2), if wb >= 0.5 { rg } else { 0 }, ta));
            }
            let (dz, rg) = self.phase_f(&bc[..nb], h2, ht, ws);
            return Some((h2 + dz, if wb >= 0.5 { rg } else { 0 }, ta));
        }
        let (mut dsum, mut wb, mut best) = (0.0, 0.0, (0.0f64, 0u8));
        let mut tsum = 0.0;
        let (mut ssum, mut kref) = (0.0, usize::MAX);
        for &(gx, gy, wgt) in &cs {
            let k = self.corner[self.clampi(gx, gy)];
            if k == u32::MAX || wgt <= 0.0 { continue; }
            let ck = &self.corners[k as usize];
            let (fz, rg) = ck.f(self.step_m, h2);
            let g = if ws > 0.0 { self.edge_g(ws, ck.tread_e(self.step_m, h2), fz - h2) } else { 1.0 };
            dsum += g * wgt * (fz - h2);
            wb += wgt;
            tsum += wgt * self.corner_ta[k as usize];
            if self.ends_label { ssum += g * wgt * (ck.f_slope(self.step_m, h2) - 1.0); }
            if wgt > best.0 { best = (wgt, rg); kref = k as usize; }
        }
        if self.ends_label && wb > 0.0 && kref != usize::MAX { best.1 = limb_of_slope(1.0 + ssum / wb, self.corners[kref].limb_slopes(), best.1); }
        Some((h2 + dsum, if wb >= 0.5 { best.1 } else { 0 }, if wb > 0.0 { tsum / wb } else { 0.0 }))
    }

    /// Lane PHASE (SHARP_BLEND_PHASE_RULE): the height change at a 2 m post, sum over the bare corners `bc` (corner index,
    /// bilinear weight) of w (Fk*(h2) - h2), and the limb of the heaviest corner. Fk* is corner k's F on its blended ladder: each
    /// level L of k's ladder becomes the weighted mean of L and the other corners' nearest levels n_j, each weighted w_j x
    /// agree(|n_j - L|) (1 within the merge tolerance, fading linearly to 0 at twice it); the knot fractions are blended the same way
    /// (agreement at both boundaries of the period). Corners that agree get the same Fk* (one F per post, no limb mixing); corners that
    /// disagree keep their own F (the F-value blend); the result is continuous in position (no reference corner).
    /// Lane TREADSOIL: `ws` = the cell's soil weight at the post (0 when ct_sharp_soil_edge is off): each corner's term takes the
    /// factor edge_g(ws, e_k), e_k from the corner's blended knots (2-limb corners: 0).
    /// Lane CELL-RISER: `ht` = the height the corners' F are evaluated at (h2 + dh; = h2 when ct_sharp_cell_riser is off).
    fn phase_f(&self, bc: &[(usize, f64)], h2: f64, ht: f64, ws: f64) -> (f64, u8) {
        let step = self.step_m;
        let tol = self.tol;
        let agree = |d: f64| if d <= tol { 1.0 } else if d >= 2.0 * tol { 0.0 } else { (2.0 * tol - d) / tol };
        let near = |c: &Corner, v: f64| { let (a, b) = c.period(step, v); if v - a <= b - v { a } else { b } };
        let (mut dz, mut best) = (0.0, (0.0f64, 0u8));
        let (mut ssum, mut wsum, mut kref) = (0.0, 0.0, 0usize);
        for (m, &(k, wk)) in bc.iter().enumerate() {
            let ck = &self.corners[k];
            let bl = |l: f64| {
                let (mut num, mut den) = (wk * l, wk);
                for (j2, &(j, wj)) in bc.iter().enumerate() {
                    if j2 == m { continue; }
                    let nj = near(&self.corners[j], l);
                    let a = wj * agree((nj - l).abs());
                    num += a * nj; den += a;
                }
                num / den
            };
            let (mut s0, mut s1) = ck.period(step, ht);
            let mut dir = 0i8;
            let mut res: Option<(f64, u8, f64, f64)> = None;
            // lane PERIOD-ENDS-2: a massive corner is the identity outside its [s0, s1] (its own F below)
            let massive_out = ck.massive_out(ht);
            for _ in 0..if massive_out { 0 } else { 7 } {
                let (b0, b1) = (bl(s0), bl(s1));
                if ht < b0 { if dir > 0 { break; } dir = -1; (s0, s1) = ck.period(step, s0 - 1e-3); continue; }
                if ht > b1 { if dir < 0 { break; } dir = 1; (s0, s1) = ck.period(step, s1 + 1e-3); continue; }
                let (uk, fk) = ck.knots();
                let (mut u, mut f, mut den) = ([wk * uk[0], wk * uk[1]], [wk * fk[0], wk * fk[1]], wk);
                for (j2, &(j, wj)) in bc.iter().enumerate() {
                    if j2 == m { continue; }
                    let cj = &self.corners[j];
                    let a = wj * agree((near(cj, s0) - s0).abs()).min(agree((near(cj, s1) - s1).abs()));
                    if a <= 0.0 { continue; }
                    // lane SPLIT-2M ct_sharp_blend_period: j's levels near s0 and s1 must bound one period of j's ladder
                    if self.same_period && cj.period(step, 0.5 * (s0 + s1)) != (near(cj, s0), near(cj, s1)) { continue; }
                    let (uj, fj) = cj.knots();
                    u[0] += a * uj[0]; u[1] += a * uj[1]; f[0] += a * fj[0]; f[1] += a * fj[1]; den += a;
                }
                let t = b1 - b0;
                let hk = [b0, b0 + t * u[0] / den, b0 + t * u[1] / den, b1];
                let zk = [b0, b0 + t * f[0] / den, b0 + t * f[1] / den, b1];
                let seg = |a: usize| { let dh = hk[a + 1] - hk[a]; if dh > 1e-9 { zk[a] + (ht - hk[a]) * (zk[a + 1] - zk[a]) / dh } else { zk[a + 1] } };
                let e = if ws > 0.0 && ck.has_tread() { edge_e(ht - hk[0], hk[1] - hk[0], hk[2] - hk[0]) } else { 0.0 };
                let sl = |a: usize| { let dh = hk[a + 1] - hk[a]; if dh > 1e-9 { (zk[a + 1] - zk[a]) / dh } else { 1.0 } };
                res = Some(if ht < hk[1] { (seg(0), 3, e, sl(0)) } else if ht < hk[2] { (seg(1), 1, e, sl(1)) } else { (seg(2), 2, e, sl(2)) });
                break;
            }
            // the blended period walk did not converge: this corner's own F
            let (fz, rg, e, sl) = res.unwrap_or_else(|| { let (fz, rg) = ck.f(step, ht); (fz, rg, if ws > 0.0 { ck.tread_e(step, ht) } else { 0.0 }, if self.ends_label { ck.f_slope(step, ht) } else { 1.0 }) });
            let g = self.edge_g(ws, e, fz - h2);
            dz += g * wk * (fz - h2);
            ssum += g * wk * (sl - 1.0);
            wsum += wk;
            if wk > best.0 { best = (wk, rg); kref = k; }
        }
        if self.ends_label && wsum > 0.0 { best.1 = limb_of_slope(1.0 + ssum / wsum, self.corners[kref].limb_slopes(), best.1); }
        (dz, best.1)
    }
}

/// Lane SHARP2: the longest period (x the step period) a 3-limb ladder accepts: pf-region var 2's default unit_max (1500
/// permille: packages up to 1.5 P are one period, the closure subdivides longer ones).
pub const TREAD_PERIOD_MAX: f64 = 1.5;

/// Lane SHARP2: the ladder (ascending period boundaries) of a 3-limb T1b post. `cand`: the bare neighbours' boundary estimates
/// (exact ones only) sorted ascending, the post's own (tread / cliff posts) flagged; `own` its s0 (tread, regime 3) or s1 (cliff,
/// regime 2), NaN for a talus post (regime 1), whose period is solved here. None: no consistent talus period.
#[allow(clippy::too_many_arguments)]
fn ladder3(cand: &[(f64, bool)], h: f64, d: f64, own: f64, rg: u8, r: f64, t3: Tread3, step: f64, tol: f64, reach: f64, consensus: bool, st: &mut [u64; 7]) -> Option<Vec<f64>> {
    let mut lad: Vec<f64> = Vec::new();
    let mut own = own;
    let mut s = 0usize;
    while s < cand.len() {
        let mut e = s + 1;
        while e < cand.len() && cand[e].0 - cand[e - 1].0 < tol { e += 1; }
        let m = e - s;
        let med = if m % 2 == 1 { cand[s + m / 2].0 } else { 0.5 * (cand[s + m / 2 - 1].0 + cand[s + m / 2].0) };
        if !cand[s..e].iter().any(|v| v.1) {
            lad.push(med);
        } else if consensus && m > 1 && ((rg == 3 && med <= h) || (rg == 2 && med >= h)) {
            // ct_sharp_consensus: the own boundary takes its cluster's median (neighbouring corners share the level)
            own = med;
        }
        s = e;
    }
    let (ut, ub) = t3.breaks(r);
    let _ = reach;
    let tmax = TREAD_PERIOD_MAX * step + tol;
    let srt = |l: &mut Vec<f64>| l.sort_by(|a, b| a.total_cmp(b));
    match rg {
        3 | 2 => {
            let (lo, hi) = (own.min(h), own.max(h));
            lad.retain(|v| !(lo - tol < *v && *v < hi + tol));
            lad.push(own);
            srt(&mut lad);
            // a period longer than tmax is not one period (pf-region var 2: packages up to unit_max x P are one period, longer ones
            // are subdivided; the closure keeps periods within [unit_min, unit_max] x P): such a level is not this post's boundary
            if rg == 3 {
                lad.retain(|v| !(*v > own && *v - own > tmax));
                loop {
                    let Some(&up) = lad.iter().find(|v| **v > own) else { break };
                    if h <= own + ut * (up - own) + 1e-9 { break; }
                    lad.retain(|v| *v != up); st[5] += 1;
                }
                // no boundary above within the radius and the period cap: the regular period (at least long enough to hold h on
                // the tread). Unlike ct_sharp_v1 no `reach` bound: at 32 m the cliff limb is narrower than a post and treads are
                // 1-2 posts wide, so an unseen boundary says little about its distance
                if !lad.iter().any(|v| *v > own) { lad.push(own + step.max((h - own) / ut * (1.0 + 1e-9))); st[6] += 1; }
            } else {
                lad.retain(|v| !(*v < own && own - *v > tmax));
                loop {
                    let Some(&dn) = lad.iter().rev().find(|v| **v < own) else { break };
                    if h >= dn + ub * (own - dn) - 1e-9 { break; }
                    lad.retain(|v| *v != dn); st[5] += 1;
                }
                if !lad.iter().any(|v| *v < own) { lad.insert(0, own - step.max((own - h) / (1.0 - ub).max(1e-6) * (1.0 + 1e-9))); st[6] += 1; }
            }
        }
        _ => {
            // talus limb: z - h = (r - 1)(h - s0) - T k1
            let k1 = t3.ft * (r / t3.p - 1.0);
            let ok = |s0: f64, t: f64| t > 0.0 && t <= tmax && { let u = (h - s0) / t; u >= ut - 1e-9 && u <= ub + 1e-9 };
            let mut sol: Option<(f64, f64, usize)> = None;
            for &l in lad.iter().filter(|v| **v > h).take(2) {
                let s0 = (d - (r - 1.0) * h + l * k1) / (k1 + 1.0 - r);
                if ok(s0, l - s0) { sol = Some((s0, l - s0, 0)); break; }
            }
            if sol.is_none() {
                for &l in lad.iter().rev().filter(|v| **v < h).take(2) {
                    let t = ((r - 1.0) * (h - l) - d) / k1;
                    if ok(l, t) { sol = Some((l, t, 1)); break; }
                }
            }
            if sol.is_none() { let s0 = h - (d + step * k1) / (r - 1.0); if ok(s0, step) { sol = Some((s0, step, 2)); } }
            if sol.is_none() {
                let um = 0.5 * (ut + ub);
                let t = d / ((r - 1.0) * um - k1);
                if t > 0.0 && t.is_finite() { sol = Some((h - um * t, t, 3)); }
            }
            let (s0, t, how) = sol?;
            st[how] += 1;
            let s1 = s0 + t;
            lad.retain(|v| !(s0 - tol < *v && *v < s1 + tol));
            lad.push(s0);
            lad.push(s1);
            srt(&mut lad);
        }
    }
    Some(lad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_region::t1b_ct::{remap, CtParamsV1, Geo, CT_LAYERED};

    /// A plane h = h0 + G x over an n x n lattice, staircased by pf-region's own remap (regular period, no apron), heights and
    /// fields quantized to mm / permille as T1b publishes them.
    fn plane(n: usize, g: f64, ta_mult: f64) -> (Vec<f64>, Vec<f64>, Vec<u8>, Vec<u16>, CtParamsV1, f64, f64, f64) {
        let c = CtParamsV1::v1();
        let tb = c.classes[3].talus_tan_permille as f64 / 1000.0;
        let ta = c.classes[3].cliff_tan_permille as f64 / 1000.0 * ta_mult;
        let cf = ((1.0 / tb - 1.0 / g) / (1.0 / tb - 1.0 / ta)).clamp(0.0, 1.0);
        let geo = Geo { base: f64::NAN, crest: f64::NAN };
        let (mut h, mut z, mut rg, mut cp) = (vec![0.0; n * n], vec![0.0; n * n], vec![0u8; n * n], vec![0u16; n * n]);
        for j in 0..n { for i in 0..n {
            let k = j * n + i;
            let hh = 3.0 + g * i as f64 * crate::T1B_M;
            let (zz, r) = remap(&c, &[], CT_LAYERED, tb, ta, cf, g, geo, hh, f64::NEG_INFINITY, f64::INFINITY);
            h[k] = (hh * 1000.0).round() / 1000.0; z[k] = (zz * 1000.0).round() / 1000.0; rg[k] = r; cp[k] = (cf * 1000.0).round().clamp(1.0, 1000.0) as u16;
        } }
        (h, z, rg, cp, c, tb, ta, cf)
    }

    fn inputs<'a>(n: usize, h: &'a [f64], z: &'a [f64], rg: &'a [u8], cp: &'a [u16], lit: &'a [u8], c: &CtParamsV1) -> CtInputs<'a> {
        CtInputs { n, gx0: 0, gy0: 0, h, z, regime: rg, cliff_permille: cp, litho: lit,
            classes: c.classes.iter().map(|k| CtClass { mode: k.mode, tb: k.talus_tan_permille as f64 / 1000.0, ta: k.cliff_tan_permille as f64 / 1000.0 }).collect(),
            step_m: c.step_mm as f64 / 1000.0, rockmass: (0.6, 1.5), tread: None, ends: None }
    }

    /// Agreement with pf-region's F_c: on a uniform slope, the 2 m re-profile of the B1 bicubic of the published staircase equals
    /// pf-region's remap evaluated at 2 m (the steps applied once, at the class angles). Tolerance: the published cliff fraction is
    /// rounded to permille, so the recovered cliff tangent carries up to ~1 % error and the riser foot moves by up to ~1 % of the
    /// riser height c P; bound 1.5 % of c P + 2 cm.
    #[test]
    fn plane_matches_pf_region_remap_at_2m() {
        for &(g, tm) in &[(1.0f64, 1.0f64), (1.4, 0.8), (0.8, 1.3)] {
            let n = 48;
            let (h, z, rg, cp, c, tb, ta, cf) = plane(n, g, tm);
            let lit = vec![3u8; n * n];
            let sf = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
            assert!(sf.corners.len() > 0);
            let geo = Geo { base: f64::NAN, crest: f64::NAN };
            let mut worst = 0.0f64;
            let y = 24.0 * crate::T1B_M;
            // interior 2 m posts (away from the lattice edge, where the bicubic stencil clamps)
            for k in (8 * 16)..(40 * 16) {
                let x = k as f64 * crate::POST_M;
                // the pipeline height with no 2 m detail: the B1 bicubic of the published staircase
                let fxx = x / crate::T1B_M;
                let (i0, t) = (fxx.floor() as i64, fxx - fxx.floor());
                let w = cr(t);
                let zb: f64 = (0..4).map(|m| w[m] * z[24 * n + (i0 - 1 + m as i64) as usize]).sum();
                let (zn, _) = sf.apply(x, y, zb).unwrap();
                let (want, _) = remap(&c, &[], CT_LAYERED, tb, ta, cf, g, geo, 3.0 + g * x, f64::NEG_INFINITY, f64::INFINITY);
                worst = worst.max((zn - want).abs());
            }
            let bound = 0.015 * cf * c.step_mm as f64 / 1000.0 + 0.02;
            assert!(worst < bound, "g {g} ta x{tm}: worst |2 m - F_c| {worst} (bound {bound})");
            eprintln!("g {g} ta x{tm}: c {cf:.4} worst |2 m - F_c| {worst:.4} m (bound {bound:.3})");
        }
    }

    /// The re-profile reproduces T1b's height at every T1b post (bilinear weight 1 on the post's own ladder).
    #[test]
    fn t1b_posts_reproduced() {
        let n = 32;
        let (h, z, rg, cp, c, ..) = plane(n, 1.1, 1.0);
        let lit = vec![3u8; n * n];
        let sf = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        for j in 4..n - 4 { for i in 4..n - 4 {
            let k = j * n + i;
            let (x, y) = (i as f64 * crate::T1B_M, j as f64 * crate::T1B_M);
            let (zn, _) = sf.apply(x, y, z[k]).unwrap();
            assert!((zn - z[k]).abs() < 1e-6, "post {i},{j}: {zn} vs {}", z[k]);
        } }
    }

    /// 2 m slopes: risers at the cliff tangent and benches at the talus tangent (not the smeared 32 m angles).
    #[test]
    fn two_metre_slopes_are_bimodal_at_the_class_angles() {
        let n = 48;
        let (h, z, rg, cp, c, tb, ta, _) = plane(n, 1.0, 1.0);
        let lit = vec![3u8; n * n];
        let sf = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        let y = 24.0 * crate::T1B_M;
        let zat = |x: f64| { let fxx = x / crate::T1B_M; let (i0, t) = (fxx.floor() as i64, fxx - fxx.floor()); let w = cr(t);
            let zb: f64 = (0..4).map(|m| w[m] * z[24 * n + (i0 - 1 + m as i64) as usize]).sum(); sf.apply(x, y, zb).unwrap().0 };
        let (mut n_tb, mut n_ta, mut tot) = (0, 0, 0);
        for k in (8 * 16)..(40 * 16) {
            let x = k as f64 * crate::POST_M;
            let s = (zat(x + crate::POST_M) - zat(x)) / crate::POST_M;
            tot += 1;
            if (s - tb).abs() < 0.02 { n_tb += 1; }
            if (s - ta).abs() < 0.1 { n_ta += 1; }
        }
        assert!(n_tb + n_ta >= tot * 9 / 10, "talus {n_tb} cliff {n_ta} of {tot}");
        assert!(n_ta > 0 && n_tb > 0);
    }

    /// Soil (regime 0) everywhere: untouched; F continuous in height across period boundaries.
    #[test]
    fn soil_untouched_and_f_continuous() {
        let n = 16;
        let h: Vec<f64> = (0..n * n).map(|k| (k % n) as f64 * 10.0).collect();
        let z = h.clone();
        let (rg, cp, lit) = (vec![0u8; n * n], vec![0u16; n * n], vec![3u8; n * n]);
        let c = CtParamsV1::v1();
        let sf = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        assert!(sf.apply(100.0, 100.0, 37.25).is_none());
        let cn = Corner { r: 0.7, q: 5.0, ladder: vec![10.0, 100.0, 250.0], t3: None, ex: None };
        let mut prev = cn.f(128.0, -400.0).0;
        let mut x = -400.0;
        while x < 800.0 { x += 0.01; let v = cn.f(128.0, x).0; assert!(v >= prev - 1e-9 && v - prev < 0.06, "F jump at {x}: {prev} -> {v}"); prev = v; }
        for b in [10.0, 100.0, 250.0] { assert!((cn.f(128.0, b).0 - b).abs() < 1e-9); }
    }

    #[test]
    fn absent_key_not_serialized() {
        let p = crate::DetailParamsV1::v3();
        let v = serde_json::to_value(&p).unwrap();
        assert!(v.get("ct_sharp").is_none());
        let mut q = crate::DetailParamsV1::v3();
        q.set("ct_sharp", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap().get("ct_sharp").is_some());
    }

    /// Lane TREAD: ct_sharp_stage needs ct_sharp; post (the default) is not serialized; pre is.
    #[test]
    fn stage_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_stage", "pre").is_err());
        q.set("ct_sharp", "1").unwrap();
        let v0 = serde_json::to_value(&q).unwrap();
        q.set("ct_sharp_stage", "post").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap(), v0);
        assert!(v0["ct_sharp"].get("stage").is_none());
        q.set("ct_sharp_stage", "pre").unwrap();
        assert!(q.ct_sharp.as_ref().unwrap().pre());
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("stage").is_some());
        assert!(q.set("ct_sharp_stage", "mid").is_err());
    }

    /// Lane SHARP2: a plane staircased by pf-region's 3-limb remap3 (ct_tread_v1, closure eta, regular period, local gradient S =
    /// g), published as T1b does: heights mm, ct_cliff_permille = fc, ct_tread_permille = ft, ct_tread_tan_permille = tt.
    #[allow(clippy::type_complexity)]
    fn plane3(n: usize, g: f64) -> (Vec<f64>, Vec<f64>, Vec<u8>, Vec<u16>, Vec<u16>, Vec<u16>, CtParamsV1, Box<dyn Fn(f64) -> f64>) {
        use pf_region::t1b_ct::{remap3, split, CtTreadParamsV1, CtV2ParamsV1};
        let mut c = CtParamsV1::v1();
        c.apply = 2;
        let mut v = CtV2ParamsV1::v1();
        v.tread = Some(CtTreadParamsV1::v1());
        c.v2 = Some(v);
        let t = c.v2.as_ref().unwrap().tread.clone().unwrap();
        let tt = t.tread_tan_permille as f64 / 1000.0;
        let (tb, ta, cf) = split(&c, 3, g, 1.0).unwrap();
        let geo = Geo { base: f64::NAN, crest: f64::NAN };
        let (mut h, mut z, mut rg, mut cp, mut ftp, mut ttp) = (vec![0.0; n * n], vec![0.0; n * n], vec![0u8; n * n], vec![0u16; n * n], vec![0u16; n * n], vec![0u16; n * n]);
        let c1 = c.clone();
        let f = move |hh: f64| remap3(&c1, &t, &[], CT_LAYERED, tb, ta, cf, g, geo, hh, f64::NEG_INFINITY, f64::INFINITY, tt, f64::NAN, &mut |_, _| None);
        for j in 0..n { for i in 0..n {
            let k = j * n + i;
            let hh = 3.0 + g * i as f64 * crate::T1B_M;
            let r = f(hh);
            assert!(r.limbs3);
            h[k] = (hh * 1000.0).round() / 1000.0; z[k] = (r.z * 1000.0).round() / 1000.0; rg[k] = r.rg;
            cp[k] = (r.fc * 1000.0).round().clamp(1.0, 1000.0) as u16; ftp[k] = (r.ft * 1000.0).round() as u16; ttp[k] = (tt * 1000.0).round() as u16;
        } }
        (h, z, rg, cp, ftp, ttp, c, Box::new(move |hh| f(hh).z))
    }

    fn bicubic_row(z: &[f64], n: usize, row: usize, x: f64) -> f64 {
        let fxx = x / crate::T1B_M;
        let (i0, t) = (fxx.floor() as i64, fxx - fxx.floor());
        let w = cr(t);
        (0..4).map(|m| w[m] * z[row * n + (i0 - 1 + m as i64) as usize]).sum()
    }

    /// Lane SHARP2: agreement with pf-region's 3-limb remap3 at 2 m on uniform slopes (tread, talus and cliff applied once, at
    /// the published angles), and every T1b post reproduced. Tolerance: ft, fc and tt are published to permille, so the
    /// recovered cliff tangent and the limb breakpoints carry permille-level errors; the error is vertical on limbs up to 45 deg
    /// and horizontal (vertical / slope) on steeper ones; bound 0.5 % of the period + 5 cm.
    #[test]
    fn tread_plane_matches_pf_region_remap3_at_2m() {
        for &g in &[0.85f64, 1.0, 1.3] {
            let n = 48;
            let (h, z, rg, cp, ftp, ttp, c, want) = plane3(n, g);
            let lit = vec![3u8; n * n];
            let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
            ci.tread = Some((&ftp, &ttp));
            let sf = SharpField::build(&ci, &SharpParams::default());
            assert!(sf.tread && sf.corners.iter().all(|k| k.t3.is_some()) && sf.corners.len() > 1000, "g {g}: {}", sf.stats);
            // the cliff limb is narrower than a 32 m post: a row may hold no T1b cliff post (the 2 m test below sees it)
            assert!(rg.iter().any(|&v| v == 3) && rg.iter().any(|&v| v == 1), "g {g}");
            let y = 24.0 * crate::T1B_M;
            let mut worst = 0.0f64;
            let mut wx = (0.0, 0u8, 0.0);
            for k in (8 * 16)..(40 * 16) {
                let x = k as f64 * crate::POST_M;
                let (zn, rr) = sf.apply(x, y, bicubic_row(&z, n, 24, x)).unwrap();
                // vertical error, as a horizontal offset on limbs steeper than 45 deg (a riser foot moved by 0.2 m is 1.1 m of
                // height on an 80 deg face)
                let sl = ((want(3.0 + g * (x + 0.5)) - want(3.0 + g * (x - 0.5))).abs()).max(1.0);
                let e = (zn - want(3.0 + g * x)).abs() / sl;
                if e > worst { worst = e; wx = (x, rr, zn); }
            }
            let bound = 0.005 * c.step_mm as f64 / 1000.0 + 0.05;
            eprintln!("g {g}: worst |2 m - remap3| (vertical / max(1, slope)) {worst:.4} m at x {:.0} m (2 m regime {}, z {:.2}) (bound {bound:.3}); T1b posts 3-limb {}, talus periods above/below/step/centred {}/{}/{}/{}", wx.0, wx.1, wx.2,
                sf.stats["tread"]["posts_3limb"], sf.stats["tread"]["talus_period_from_level_above"], sf.stats["tread"]["talus_period_from_level_below"], sf.stats["tread"]["talus_period_step"], sf.stats["tread"]["talus_period_centred"]);
            assert!(worst < bound, "g {g}: worst {worst} (bound {bound})");
            for j in 4..n - 4 { for i in 4..n - 4 {
                let k = j * n + i;
                let (zn, _) = sf.apply(i as f64 * crate::T1B_M, j as f64 * crate::T1B_M, z[k]).unwrap();
                assert!((zn - z[k]).abs() < 1e-6, "g {g} post {i},{j} (regime {}): {zn} vs {}", rg[k], z[k]);
            } }
        }
    }

    /// Lane SHARP2: 2 m slopes fall on the three published angles (tread tt, talus tb, cliff ta) and the 2 m regime marks them
    /// (3 tread, 1 talus, 2 cliff); the tread is below 20 deg.
    #[test]
    fn tread_two_metre_slopes_are_trimodal() {
        let n = 48;
        let g = 1.0;
        let (h, z, rg, cp, ftp, ttp, c, _) = plane3(n, g);
        let lit = vec![3u8; n * n];
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.tread = Some((&ftp, &ttp));
        let sf = SharpField::build(&ci, &SharpParams::default());
        let tt = ttp[0] as f64 / 1000.0;
        let tb = c.classes[3].talus_tan_permille as f64 / 1000.0;
        assert!(tt < 20f64.to_radians().tan());
        let y = 24.0 * crate::T1B_M;
        let at = |x: f64| sf.apply(x, y, bicubic_row(&z, n, 24, x)).unwrap();
        let (mut n3, mut n1, mut n2, mut tot) = (0, 0, 0, 0);
        for k in (8 * 16)..(40 * 16) {
            let x = k as f64 * crate::POST_M;
            let (a, b) = (at(x), at(x + crate::POST_M));
            let s = (b.0 - a.0) / crate::POST_M;
            tot += 1;
            if (s - tt).abs() < 0.01 && a.1 == 3 { n3 += 1; }
            if (s - tb).abs() < 0.02 && a.1 == 1 { n1 += 1; }
            if s > 3.0 && a.1 == 2 { n2 += 1; }
        }
        eprintln!("tread {n3} talus {n1} cliff {n2} of {tot}");
        assert!(n3 > 0 && n1 > 0 && n2 > 0);
        assert!(n3 + n1 + n2 >= tot * 85 / 100, "tread {n3} talus {n1} cliff {n2} of {tot}");
    }

    /// Lane SHARP2: tread fields present but zero everywhere (no 3-limb post) give exactly the 2-limb re-profile.
    #[test]
    fn zero_tread_fields_are_the_two_limb_profile() {
        let n = 32;
        let (h, z, rg, cp, c, ..) = plane(n, 1.1, 1.0);
        let lit = vec![3u8; n * n];
        let zero = vec![0u16; n * n];
        let a = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.tread = Some((&zero, &zero));
        let b = SharpField::build(&ci, &SharpParams::default());
        assert!(!a.tread && b.tread && a.stats.get("tread").is_none());
        for k in 0..400 {
            let (x, y) = (150.0 + 1.3 * k as f64, 200.0 + 0.7 * k as f64);
            let zz = 3.0 + 1.1 * x - 4.0;
            assert_eq!(a.apply_ta(x, y, zz), b.apply_ta(x, y, zz));
        }
    }

    /// Lane SHARP2: ct_sharp_consensus needs ct_sharp, is not serialized when off; with it the 3-limb plane still matches remap3.
    #[test]
    fn consensus_key_opt_in_and_plane() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_consensus", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        let v0 = serde_json::to_value(&q).unwrap();
        assert!(v0["ct_sharp"].get("t3_consensus").is_none());
        q.set("ct_sharp_consensus", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"]["t3_consensus"] == serde_json::json!(true));
        let n = 48;
        let (h, z, rg, cp, ftp, ttp, c, want) = plane3(n, 1.0);
        let lit = vec![3u8; n * n];
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.tread = Some((&ftp, &ttp));
        let sp = SharpParams { t3_consensus: true, ..SharpParams::default() };
        let sf = SharpField::build(&ci, &sp);
        let y = 24.0 * crate::T1B_M;
        let mut worst = 0.0f64;
        for k in (8 * 16)..(40 * 16) {
            let x = k as f64 * crate::POST_M;
            let (zn, _) = sf.apply(x, y, bicubic_row(&z, n, 24, x)).unwrap();
            let sl = ((want(3.0 + (x + 0.5)) - want(3.0 + (x - 0.5))).abs()).max(1.0);
            worst = worst.max((zn - want(3.0 + x)).abs() / sl);
        }
        assert!(worst < 0.005 * c.step_mm as f64 / 1000.0 + 0.05, "consensus worst {worst}");
    }

    /// Lane TREAD: apply_ta returns apply's height and regime, and the cliff tangent of the bare corners.
    #[test]
    fn apply_ta_matches_apply() {
        let n = 24;
        let (h, z, rg, cp, c, _tb, ta, _cf) = plane(n, 1.1, 1.0);
        let lit = vec![3u8; n * n];
        let sf = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        let mut seen = 0;
        for k in 0..200 {
            let (x, y) = (200.0 + 1.7 * k as f64, 300.0 + 0.9 * k as f64);
            let zz = 3.0 + 1.1 * x - 5.0;
            match (sf.apply(x, y, zz), sf.apply_ta(x, y, zz)) {
                (None, None) => {}
                (Some(a), Some(b)) => { assert!(a.0 == b.0 && a.1 == b.1); assert!((b.2 / ta - 1.0).abs() < 0.02, "ta {} vs {ta}", b.2); seen += 1; }
                _ => panic!("apply / apply_ta disagree"),
            }
        }
        assert!(seen > 100);
    }

    fn phase_params() -> SharpParams { SharpParams { blend: Some(SHARP_BLEND_PHASE_RULE.to_string()), ..SharpParams::default() } }

    /// Lane PHASE: ct_sharp_blend needs ct_sharp, is not serialized when off (or =f), accepts phase / nearest only.
    #[test]
    fn blend_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_blend", "phase").is_err());
        q.set("ct_sharp", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("blend").is_none());
        q.set("ct_sharp_blend", "f").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("blend").is_none());
        assert!(q.set("ct_sharp_blend", "bilinear").is_err());
        q.set("ct_sharp_blend", "phase").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["blend"], serde_json::json!(SHARP_BLEND_PHASE_RULE));
        q.set("ct_sharp_blend", "nearest").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["blend"], serde_json::json!(SHARP_BLEND_NEAREST_RULE));
    }

    /// Lane PHASE: on pf-region's 3-limb plane (consistent corners) the phase blend still matches remap3 at 2 m within the
    /// SHARP2 bound and reproduces every T1b post exactly; on the 2-limb plane it matches the F-value blend.
    #[test]
    fn phase_plane_matches_remap3_and_posts() {
        for &g in &[0.85f64, 1.0, 1.3] {
            let n = 48;
            let (h, z, rg, cp, ftp, ttp, c, want) = plane3(n, g);
            let lit = vec![3u8; n * n];
            let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
            ci.tread = Some((&ftp, &ttp));
            let sf = SharpField::build(&ci, &phase_params());
            assert_eq!(sf.blend, 1);
            let y = 24.0 * crate::T1B_M;
            let mut worst = 0.0f64;
            for k in (8 * 16)..(40 * 16) {
                let x = k as f64 * crate::POST_M;
                let (zn, _) = sf.apply(x, y, bicubic_row(&z, n, 24, x)).unwrap();
                let sl = ((want(3.0 + g * (x + 0.5)) - want(3.0 + g * (x - 0.5))).abs()).max(1.0);
                worst = worst.max((zn - want(3.0 + g * x)).abs() / sl);
            }
            let bound = 0.005 * c.step_mm as f64 / 1000.0 + 0.05;
            eprintln!("phase g {g}: worst {worst:.4} (bound {bound:.3})");
            assert!(worst < bound, "g {g}: worst {worst}");
            for j in 4..n - 4 { for i in 4..n - 4 {
                let k = j * n + i;
                let (zn, _) = sf.apply(i as f64 * crate::T1B_M, j as f64 * crate::T1B_M, z[k]).unwrap();
                assert!((zn - z[k]).abs() < 1e-6, "g {g} post {i},{j}: {zn} vs {}", z[k]);
            } }
        }
        let n = 24;
        let (h, z, rg, cp, c, ..) = plane(n, 1.1, 1.0);
        let lit = vec![3u8; n * n];
        let a = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        let b = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &phase_params());
        let mut worst = 0.0f64;
        for k in 0..300 {
            let (x, y) = (200.0 + 1.3 * k as f64, 250.0 + 0.7 * k as f64);
            let zz = 3.0 + 1.1 * x - 4.0;
            if let (Some(p), Some(q)) = (a.apply(x, y, zz), b.apply(x, y, zz)) { worst = worst.max((p.0 - q.0).abs()); }
        }
        eprintln!("2-limb plane: worst |phase - F blend| {worst:.4} m");
        assert!(worst < 0.5, "2-limb plane phase vs F blend {worst}");
    }

    /// Lane PHASE: two corner columns whose ladders are offset by dl (the SHARP2 mixed-cell case). F-value blend: inside the cell
    /// one corner is already on its talus while the other is still on its tread, so posts labelled tread climb at talus-like
    /// slopes. Phase blend with dl within the merge tolerance: one F; the tread limb climbs at p G + (1 - p) dl / 32 at most. For
    /// every dl (agreeing, fading, disagreeing ladders): Fb <= h, monotone in h, and continuous along the row across the cell
    /// mid-line (x 48 m), where a reference-corner switch would put a seam.
    #[test]
    fn phase_blend_keeps_offset_benches_flat() {
        let n = 4;
        let t3 = Tread3 { p: 0.3, ft: 0.1, fb: 0.5, fc: 0.4 };
        let r = 0.9;
        let ub = t3.breaks(r).1;
        let q = t3.fc / (1.0 - ub);
        for &dl in &[5.0f64, 9.0, 20.0] {
            let lad = |o: f64| (0..6).map(|k| o + 100.0 * k as f64).collect::<Vec<f64>>();
            let corners = vec![Corner { r, q, ladder: lad(0.0), t3: Some(t3), ex: None }, Corner { r, q, ladder: lad(dl), t3: Some(t3), ex: None }];
            let corner: Vec<u32> = (0..n * n).map(|i| if i % n < 2 { 0 } else { 1 }).collect();
            let mk = |blend: u8| SharpField { same_period: false, n, gx0: 0, gy0: 0, step_m: 100.0, d: vec![0.0; n * n], corner: corner.clone(), corners: corners.clone(), corner_ta: vec![1.5, 1.5], tread: true, stats: serde_json::Value::Null, blend, tol: 6.0, soil_edge: 0.0, soil_tw: 0.0, soil: Vec::new(), ends_label: false, dh_post: Vec::new(), dh_mode: 0 };
            let (fv, ph) = (mk(0), mk(1));
            let y = 48.0;
            let lim = t3.p + (1.0 - t3.p) * dl / crate::T1B_M + 1e-6;
            let (mut fv_max, mut ph_max, mut nph, mut jump) = (0.0f64, 0.0f64, 0, 0.0f64);
            for &h0 in &[80.0f64, 90.0, 100.0, 110.0, 120.0] {
                let mut prev: Option<(f64, f64, u8, u8)> = None;
                for k in 0..=640 {
                    let x = 32.0 + 0.05 * k as f64;
                    let h2 = h0 + (x - 32.0);
                    let (za, ra) = fv.apply(x, y, h2).unwrap();
                    let (zb, rb) = ph.apply(x, y, h2).unwrap();
                    assert!(zb <= h2 + 1e-9, "dl {dl}: Fb > h at x {x}");
                    assert!(ph.apply(x, y, h2 + 0.01).unwrap().0 >= zb - 1e-9, "dl {dl}: not monotone at x {x}");
                    if let Some((pa, pb, pra, prb)) = prev {
                        let (sa, sb) = ((za - pa) / 0.05, (zb - pb) / 0.05);
                        jump = jump.max(sb.abs());
                        if ra == 3 && pra == 3 { fv_max = fv_max.max(sa); }
                        if rb == 3 && prb == 3 { ph_max = ph_max.max(sb); nph += 1; }
                    }
                    prev = Some((za, zb, ra, rb));
                }
            }
            eprintln!("offset {dl} m: tread slope max F-blend {fv_max:.3}, phase {ph_max:.3} (limit {lim:.3}; {nph} phase tread steps); steepest phase step {jump:.2} (cliff q {q:.2})");
            // continuity: no step steeper than the cliff limb (+ the blend's own tilt)
            assert!(jump < q * 1.5 + 1.0, "dl {dl}: discontinuity, slope {jump}");
            if dl < 6.0 {
                assert!(nph > 200);
                assert!(ph_max <= lim, "phase tread slope {ph_max} > {lim}");
                assert!(fv_max > lim, "F-value blend tread slope {fv_max} (expected the mixed-cell tilt)");
            }
        }
    }

    /// Lane TREADSOIL: ct_sharp_soil_edge needs ct_sharp, is not serialized when off (or =0); the field built without it is off.
    #[test]
    fn soil_edge_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_soil_edge", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("soil_edge").is_none());
        q.set("ct_sharp_soil_edge", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("soil_edge").is_none());
        q.set("ct_sharp_soil_edge", "1").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["soil_edge"], serde_json::json!(SHARP_SOIL_EDGE_RULE));
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("soil_edge_tw").is_none());
        assert!(q.set("ct_sharp_soil_edge_tw", "-1").is_err());
        q.set("ct_sharp_soil_edge_tw", "0.5").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["soil_edge_tw"], serde_json::json!(0.5));
        let mut q2 = crate::DetailParamsV1::v3();
        q2.set("ct_sharp", "1").unwrap();
        assert!(q2.set("ct_sharp_soil_edge_tw", "0.5").is_err());
        let n = 24;
        let (h, z, rg, cp, c, ..) = plane(n, 1.1, 1.0);
        let lit = vec![3u8; n * n];
        let a = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams::default());
        assert!(a.soil_edge == 0.0 && a.soil.is_empty());
        let b = SharpField::build(&inputs(n, &h, &z, &rg, &cp, &lit, &c), &SharpParams { soil_edge: Some(SHARP_SOIL_EDGE_RULE.to_string()), ..SharpParams::default() });
        assert!(b.soil_edge == SOIL_EDGE_E0 && b.soil.len() == n * n);
    }

    /// Lane TREADSOIL: a cell with two bare 3-limb corners (west, x 32 m) and two soil corners (east, x 64 m) on a slope G = 1.
    /// Linear soil fade (key off): posts labelled tread climb at p G + (1 - p)(h2 - s0) / 32 (the cut depth faded over the cell).
    /// Soil edge (wedge limit TW 0.5 here): where the cut F - h2 is within 8 m (the thin-wedge range) tread posts climb at p G exactly (the bench kept whole
    /// up to the mantle edge, x 48 m); a cut of 16 m or more keeps the linear fade; the colluvium wedge lies
    /// between the bench F and the soil surface h2, T1b posts are reproduced (bare: F, soil: h2) and the row is continuous; both
    /// blends (F values, phase).
    #[test]
    fn soil_edge_keeps_benches_flat_to_the_mantle_edge() {
        let n = 4;
        let t3 = Tread3 { p: 0.3, ft: 0.1, fb: 0.5, fc: 0.4 };
        let r = 0.9;
        let ub = t3.breaks(r).1;
        let q = t3.fc / (1.0 - ub);
        let lad: Vec<f64> = (0..6).map(|k| 100.0 * k as f64).collect();
        let corners = vec![Corner { r, q, ladder: lad, t3: Some(t3), ex: None }];
        let corner: Vec<u32> = (0..n * n).map(|i| if i % n < 2 { 0 } else { u32::MAX }).collect();
        let soil: Vec<bool> = (0..n * n).map(|i| i % n >= 2).collect();
        for blend in [0u8, 1] {
            let mk = |e0: f64| SharpField { same_period: false, n, gx0: 0, gy0: 0, step_m: 100.0, d: vec![0.0; n * n], corner: corner.clone(), corners: corners.clone(), corner_ta: vec![1.5], tread: true, stats: serde_json::Value::Null, blend, tol: 6.0, soil_edge: e0, soil_tw: 0.5, soil: if e0 > 0.0 { soil.clone() } else { Vec::new() }, ends_label: false, dh_post: Vec::new(), dh_mode: 0 };
            let (off, on) = (mk(0.0), mk(SOIL_EDGE_E0));
            let y = 48.0;
            let (mut off_max, mut on_max, mut non, mut jump) = (0.0f64, 0.0f64, 0, 0.0f64);
            let mut nfar = 0;
            for &h0 in &[95.0f64, 100.5, 105.0, 115.0, 125.0, 150.0, 190.0] {
                let mut prev: Option<(f64, f64, u8, u8)> = None;
                for k in 0..=640 {
                    let x = 32.0 + 0.05 * k as f64;
                    let h2 = h0 + (x - 32.0);
                    // the soil T1b post (x 64 m) is in the all-soil cell east: untouched (None)
                    let (za, ra) = off.apply(x, y, h2).unwrap_or((h2, 0));
                    let (zb, rb) = on.apply(x, y, h2).unwrap_or((h2, 0));
                    let (fz, _) = corners[0].f(100.0, h2);
                    assert!(zb <= h2 + 1e-9 && zb >= fz - 1e-9, "blend {blend}: z {zb} outside [F {fz}, h2 {h2}] at x {x}");
                    assert_eq!(ra, rb, "labels unchanged");
                    if k == 0 { assert!((zb - fz).abs() < 1e-9 && (za - fz).abs() < 1e-9, "bare T1b post"); }
                    if k == 640 { assert!((zb - h2).abs() < 1e-9, "soil T1b post"); }
                    // a thick wedge (cut >= 16 m) keeps the linear fade
                    if h2 - fz >= crate::T1B_M * 0.5 { assert!((za - zb).abs() < 1e-9, "thick wedge at x {x}: {za} vs {zb}"); nfar += 1; }
                    if let Some((pa, pb, pra, prb)) = prev {
                        let (sa, sb) = ((za - pa) / 0.05, (zb - pb) / 0.05);
                        jump = jump.max(sb.abs());
                        if ra == 3 && pra == 3 { off_max = off_max.max(sa); }
                        // (test limit TW 0.5) within the thin-wedge range (cut <= 32 m x TW x (1 - e0) = 8 m) the bench is whole
                        if rb == 3 && prb == 3 && h2 - fz <= 8.0 - 1e-6 { on_max = on_max.max(sb); non += 1; }
                    }
                    prev = Some((za, zb, ra, rb));
                }
            }
            eprintln!("blend {blend}: tread slope max linear fade {off_max:.3}, soil edge {on_max:.3} (p {}; {non} steps); steepest step {jump:.2}", t3.p);
            assert!(non > 100 && nfar > 100, "{non} {nfar}");
            assert!(on_max <= t3.p + 1e-6, "soil-edge tread slope {on_max}");
            assert!(off_max > t3.p + 0.1, "linear fade tread slope {off_max} (expected the soil tilt)");
            assert!(jump < q * 1.5 + 2.0, "discontinuity, slope {jump}");
        }
    }

    /// Lane TREADSOIL: soil posts on or next to a fixed T1b post are not soil corners for the soil edge.
    #[test]
    fn soil_edge_excludes_soil_near_fixed() {
        let n = 6;
        let mut sf = SharpField { same_period: false, n, gx0: 0, gy0: 0, step_m: 100.0, d: vec![0.0; n * n], corner: vec![u32::MAX; n * n], corners: Vec::new(), corner_ta: Vec::new(), tread: true, stats: serde_json::json!({}), blend: 1, tol: 6.0, soil_edge: SOIL_EDGE_E0, soil_tw: SOIL_EDGE_TW, soil: vec![true; n * n], ends_label: false, dh_post: Vec::new(), dh_mode: 0 };
        let mut fixed = vec![0u8; n * n];
        fixed[2 * n + 2] = 1;
        sf.exclude_soil_near_fixed(&fixed);
        for y in 0..n { for x in 0..n {
            let near = (x as i64 - 2).abs() <= 1 && (y as i64 - 2).abs() <= 1;
            assert_eq!(sf.soil[y * n + x], !near, "post {x},{y}");
        } }
        assert_eq!(sf.stats["soil_edge_excluded_near_fixed_posts"], serde_json::json!(9));
    }

    /// Lane CARRY: ct_sharp_carry_strike / ct_sharp_limb_lp need ct_sharp, are not serialized when off; sub-keys need their key.
    #[test]
    fn carry_strike_and_limb_lp_keys_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_carry_strike", "1").is_err());
        assert!(q.set("ct_sharp_limb_lp", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        let base = serde_json::to_value(&q).unwrap();
        for (k, f, sk, sf) in [("ct_sharp_carry_strike", "carry_strike", "ct_sharp_carry_strike_m", "carry_strike_m"), ("ct_sharp_limb_lp", "limb_lp", "ct_sharp_limb_lp_m", "limb_lp_m")] {
            assert!(q.set(sk, "8").is_err());
            q.set(k, "0").unwrap();
            assert_eq!(serde_json::to_value(&q).unwrap(), base);
            q.set(k, "1").unwrap();
            assert!(serde_json::to_value(&q).unwrap()["ct_sharp"][f].is_string());
            assert!(q.set(sk, "1").is_err() && q.set(sk, "40").is_err());
            q.set(sk, "12").unwrap();
            assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"][sf], serde_json::json!(12.0));
            q.set(k, "0").unwrap();
            assert_eq!(serde_json::to_value(&q).unwrap(), base);
        }
    }

    /// Lane CARRY2: ct_sharp_aa / ct_sharp_riser_noflux / ct_sharp_carry_fringe need their parents, are not serialized when off;
    /// ct_sharp_aa_w needs ct_sharp_aa and is range-checked.
    #[test]
    fn carry2_keys_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        for k in ["ct_sharp_aa", "ct_sharp_riser_noflux", "ct_sharp_carry_fringe"] { assert!(q.set(k, "1").is_err()); }
        q.set("ct_sharp", "1").unwrap();
        let base = serde_json::to_value(&q).unwrap();
        assert!(q.set("ct_sharp_aa_w", "2").is_err());
        assert!(q.set("ct_sharp_carry_fringe", "1").is_err(), "fringe needs carry_strike");
        for (k, f) in [("ct_sharp_aa", "aa"), ("ct_sharp_riser_noflux", "riser_noflux")] {
            q.set(k, "0").unwrap();
            assert_eq!(serde_json::to_value(&q).unwrap(), base);
            q.set(k, "1").unwrap();
            assert!(serde_json::to_value(&q).unwrap()["ct_sharp"][f].is_string());
            q.set(k, "0").unwrap();
            assert_eq!(serde_json::to_value(&q).unwrap(), base);
        }
        q.set("ct_sharp_aa", "1").unwrap();
        assert!(q.set("ct_sharp_aa_w", "0.2").is_err() && q.set("ct_sharp_aa_w", "3").is_err());
        q.set("ct_sharp_aa_w", "2").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["aa_w"], serde_json::json!(2.0));
        q.set("ct_sharp_aa", "0").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap(), base);
        q.set("ct_sharp_carry_strike", "1").unwrap();
        let cs = serde_json::to_value(&q).unwrap();
        q.set("ct_sharp_carry_fringe", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"]["carry_fringe"].is_string());
        q.set("ct_sharp_carry_fringe", "0").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap(), cs);
    }

    /// Lane STEEP-TALUS: ct_sharp_talus_rock needs ct_sharp, is not serialized when off (or =0); the threshold sub-key needs it.
    #[test]
    fn talus_rock_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_talus_rock", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("talus_rock").is_none());
        assert!(q.set("ct_sharp_talus_rock_deg", "50").is_err());
        q.set("ct_sharp_talus_rock", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("talus_rock").is_none());
        q.set("ct_sharp_talus_rock", "1").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["talus_rock"], serde_json::json!(SHARP_TALUS_ROCK_RULE));
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("talus_rock_deg").is_none());
        assert!(q.set("ct_sharp_talus_rock_deg", "20").is_err());
        assert!(q.set("ct_sharp_talus_rock_deg", "85").is_err());
        q.set("ct_sharp_talus_rock_deg", "50").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["talus_rock_deg"], serde_json::json!(50.0));
        q.set("ct_sharp_talus_rock", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("talus_rock_deg").is_none());
    }

    /// Lane STEEP-TALUS: a talus-labelled ramp steeper than 45 deg over 8 m on both sides becomes rock; talus at its repose, the
    /// one-post talus fringe at a cliff foot and crest (steep on one side only), posts near the ring and other regimes stay.
    #[test]
    fn talus_rock_relabels_steep_ramps_not_the_cliff_fringe() {
        // profile along x (rows identical): talus 0.75 for 20 posts, cliff 3.0 for 3 posts (regime 2), talus 0.75 for 12 posts,
        // a talus-labelled ramp 1.6 (58 deg) for 12 posts, talus 0.75 to the end; heights rise with x
        let n = 64;
        let mut prof = vec![0.0f64; n];
        let mut lab = vec![1u8; n];
        for i in 1..n {
            let t = if (20..23).contains(&i) { 3.0 } else if (35..47).contains(&i) { 1.6 } else { 0.75 };
            prof[i] = prof[i - 1] + t * crate::POST_M;
            if (20..23).contains(&i) { lab[i] = 2; }
        }
        let z: Vec<f64> = (0..n * n).map(|k| prof[k % n]).collect();
        let mut rg: Vec<u8> = (0..n * n).map(|k| lab[k % n]).collect();
        rg[30 * n + 40] = 3; // a tread-labelled post on the ramp is not relabelled
        let before = rg.clone();
        let tan45 = 45f64.to_radians().tan();
        let rl = talus_rock_relabel(n, &z, &mut rg, tan45);
        let row = 30;
        let got: Vec<usize> = (0..n).filter(|&i| rg[row * n + i] == 2 && before[row * n + i] == 1).collect();
        // ramp steps are i = 35..46 (z[i] - z[i-1] = 1.6 x 2 m); (z[i] - z[i-4]) / 8 > 1 needs >= 2 ramp steps among i-3..i
        // (i >= 36) and (z[i+4] - z[i]) / 8 > 1 needs >= 2 among i+1..i+4 (i <= 44)
        assert_eq!(got, (36..=44).filter(|&i| i != 40).collect::<Vec<_>>(), "relabelled {got:?}");
        assert_eq!(rg[30 * n + 40], 3);
        // the cliff fringe (talus posts 17..26 next to the cliff 20..22) stays talus, as does gentle talus
        for i in (0..36).chain(45..n) { if before[row * n + i] == 1 { assert_eq!(rg[row * n + i], 1, "post {i}"); } }
        // ring band (4 posts) untested; every relabelled post was talus
        for &k in &rl { assert_eq!(before[k], 1); let (i, j) = (k % n, k / n); assert!(i >= 4 && j >= 4 && i < n - 4 && j < n - 4); }
        // a 60 deg threshold (tan 1.73 > the ramp's 1.6) keeps the 58 deg ramp talus
        let mut rg2 = before.clone();
        assert!(talus_rock_relabel(n, &z, &mut rg2, 60f64.to_radians().tan()).is_empty());
    }

    /// Lane PERIOD-ENDS-2: ct_sharp_ends needs ct_sharp, is not serialized when off (or =0); the label sub-key needs it and is
    /// dropped with it.
    #[test]
    fn ends_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_ends", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("ends").is_none());
        assert!(q.set("ct_sharp_ends_label", "1").is_err());
        q.set("ct_sharp_ends", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("ends").is_none());
        q.set("ct_sharp_ends", "1").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["ends"], serde_json::json!(SHARP_ENDS_RULE));
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("ends_label").is_none());
        q.set("ct_sharp_ends_label", "1").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["ends_label"], serde_json::json!(SHARP_ENDS_LABEL_RULE));
        q.set("ct_sharp_ends", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("ends_label").is_none());
        assert_eq!(SharpParams::default().ends, None);
    }

    /// Lane PERIOD-ENDS-2 test plane: h = 3 + g x on an irregular ladder (periods 50-230 m, the var 3 / record-package case the
    /// inferred ladders fill at the regular period), F through fixed knot fractions, every T1b post publishing its own ct_ends
    /// fields (mm). Returns the window inputs, the ends and the true F.
    fn ends_plane(n: usize, g: f64) -> (Vec<f64>, Vec<f64>, Vec<u8>, Vec<u16>, Vec<u16>, Vec<u16>, crate::input::CtEnds, CtParamsV1, Box<dyn Fn(f64) -> (f64, u8)>) {
        let c = CtParamsV1::v1();
        let (u, f) = ([0.2f64, 0.7], [0.06f64, 0.5]);
        let ts = [70.0f64, 180.0, 55.0, 230.0, 95.0, 140.0, 60.0, 210.0, 120.0, 80.0, 190.0, 66.0];
        let mut lad = vec![-400.0f64];
        for k in 0..60 { lad.push(lad[k] + ts[k % ts.len()]); }
        let lad2 = lad.clone();
        let fz = move |h: f64| -> (f64, u8, usize) {
            let k = lad2.partition_point(|v| *v <= h) - 1;
            let (s0, s1) = (lad2[k], lad2[k + 1]);
            let (t, x) = (s1 - s0, (h - s0) / (s1 - s0));
            let (z, rg) = if x < u[0] { (x * f[0] / u[0], 3) } else if x < u[1] { (f[0] + (x - u[0]) * (f[1] - f[0]) / (u[1] - u[0]), 1) } else { (f[1] + (x - u[1]) * (1.0 - f[1]) / (1.0 - u[1]), 2) };
            (s0 + t * z, rg, k)
        };
        let nn = n * n;
        let (mut h, mut z, mut rg, cp, ftp, mut ttp) = (vec![0.0; nn], vec![0.0; nn], vec![0u8; nn], vec![500u16; nn], vec![60u16; nn], vec![0u16; nn]);
        let mm = |v: f64| (v * 1000.0).round() as i32;
        let mut e = crate::input::CtEnds { s0: vec![0; nn], s1: vec![0; nn], s_dn: vec![0; nn], s_up: vec![0; nn], h1: vec![0; nn], z1: vec![0; nn], h2: vec![0; nn], z2: vec![0; nn], grad_permille: vec![(g * 1000.0).round() as u16; nn], kind: vec![1; nn] };
        for j in 0..n { for i in 0..n {
            let k = j * n + i;
            let hh = 3.0 + g * i as f64 * crate::T1B_M;
            let (zz, r, p) = fz(hh);
            h[k] = (hh * 1000.0).round() / 1000.0; z[k] = (zz * 1000.0).round() / 1000.0; rg[k] = r;
            ttp[k] = ((f[0] / u[0]) * g * 1000.0).round() as u16;
            let (s0, s1) = (lad[p], lad[p + 1]);
            let t = s1 - s0;
            e.s0[k] = mm(s0); e.s1[k] = mm(s1); e.s_dn[k] = mm(lad[p - 1]); e.s_up[k] = mm(lad[p + 2]);
            e.h1[k] = mm(s0 + u[0] * t); e.z1[k] = mm(s0 + f[0] * t); e.h2[k] = mm(s0 + u[1] * t); e.z2[k] = mm(s0 + f[1] * t);
        } }
        (h, z, rg, cp, ftp, ttp, e, c, Box::new(move |hh| { let (a, b, _) = fz(hh); (a, b) }))
    }

    /// Lane PERIOD-ENDS-2: with the published ends and knots the 2 m re-profile of a plane on an IRREGULAR ladder is the true F
    /// at every 2 m post (mm knots: < 1 cm), the T1b posts are reproduced, the 2 m limbs are the true limbs, and the cliff tangent
    /// is S q; the inferred ladders (key off) miss it by metres (printed). Phase and F-value blends agree (same F at every corner).
    #[test]
    fn ends_plane_irregular_ladder_is_exact() {
        for &g in &[0.9f64, 1.2] {
            let n = 48;
            let (h, z, rg, cp, ftp, ttp, e, c, want) = ends_plane(n, g);
            let lit = vec![3u8; n * n];
            let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
            ci.tread = Some((&ftp, &ttp));
            let off = SharpField::build(&ci, &phase_params());
            ci.ends = Some(&e);
            let mut pe = phase_params();
            pe.ends = Some(SHARP_ENDS_RULE.to_string());
            let on = SharpField::build(&ci, &pe);
            let on0 = SharpField::build(&ci, &SharpParams { ends: Some(SHARP_ENDS_RULE.to_string()), ..SharpParams::default() });
            assert_eq!(on.stats["ends"]["exact_posts"], serde_json::json!(n * n));
            assert_eq!(on.stats["ends"]["held_posts_old_path"], serde_json::json!(0));
            assert!(off.stats.get("ends").is_none());
            let y = 24.0 * crate::T1B_M;
            let (mut w_on, mut w_off, mut w_lab, mut cnt) = (0.0f64, 0.0f64, 0usize, 0usize);
            for k in (8 * 16)..(40 * 16) {
                let x = k as f64 * crate::POST_M;
                let z2 = bicubic_row(&z, n, 24, x);
                let (wz, wr) = want(3.0 + g * x);
                let (za, ra) = on.apply(x, y, z2).unwrap();
                let (zb, _) = on0.apply(x, y, z2).unwrap();
                let (zo, _) = off.apply(x, y, z2).unwrap();
                w_on = w_on.max((za - wz).abs()).max((zb - wz).abs());
                w_off = w_off.max((zo - wz).abs());
                cnt += 1;
                if ra != wr { w_lab += 1; }
            }
            eprintln!("ends plane g {g}: worst |2 m - F| exact {w_on:.4} m, inferred {w_off:.2} m; limb mismatches {w_lab}/{cnt}");
            assert!(w_on < 0.01, "g {g}: exact F off by {w_on}");
            assert!(w_lab * 100 < cnt, "g {g}: {w_lab} limb mismatches");
            for j in 4..n - 4 { for i in 4..n - 4 {
                let k = j * n + i;
                let (zn, _) = on.apply(i as f64 * crate::T1B_M, j as f64 * crate::T1B_M, z[k]).unwrap();
                assert!((zn - z[k]).abs() < 0.005, "g {g} post {i},{j}: {zn} vs {}", z[k]);
            } }
            let q = (1.0 - 0.5) / (1.0 - 0.7);
            for &ta in &on.corner_ta { assert!((ta - g * q).abs() < 0.01, "ta {ta} vs {}", g * q); }
        }
    }

    /// Lane PERIOD-ENDS-2: posts flagged held / moved (kind + 128) and bare posts without a period keep the corner they get without
    /// the key (same ladder, same F); the window inputs without ct_ends leave the field as without the key.
    #[test]
    fn ends_held_posts_keep_the_inferred_corner() {
        let n = 32;
        let (h, z, rg, cp, ftp, ttp, mut e, c, _) = ends_plane(n, 1.0);
        let lit = vec![3u8; n * n];
        let held: Vec<usize> = (0..n * n).filter(|k| k % 5 == 0).collect();
        for &k in &held { e.kind[k] |= 0x80; }
        e.s0[7] = i32::MIN;
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.tread = Some((&ftp, &ttp));
        let off = SharpField::build(&ci, &phase_params());
        let mut pe = phase_params();
        pe.ends = Some(SHARP_ENDS_RULE.to_string());
        let no_fields = SharpField::build(&ci, &pe);
        assert_eq!(no_fields.corners.len(), off.corners.len());
        assert!(no_fields.stats.get("ends").is_none());
        ci.ends = Some(&e);
        let on = SharpField::build(&ci, &pe);
        assert_eq!(on.stats["ends"]["held_posts_old_path"], serde_json::json!(held.len()));
        assert_eq!(on.stats["ends"]["bare_posts_without_ends"], serde_json::json!(1));
        for &k in held.iter().chain([7usize].iter()) {
            assert_eq!(on.corner[k] == u32::MAX, off.corner[k] == u32::MAX, "post {k}");
            if off.corner[k] == u32::MAX { continue; }
            let (a, b) = (&on.corners[on.corner[k] as usize], &off.corners[off.corner[k] as usize]);
            assert!(a.ex.is_none() && a.ladder == b.ladder && a.r == b.r && a.q == b.q, "post {k}");
            assert_eq!(on.corner_ta[on.corner[k] as usize], off.corner_ta[off.corner[k] as usize]);
        }
    }

    /// Lane PERIOD-ENDS-2: apron (kind 5) and massive (kind 6) corners: the apron's own period [s0, s1] at and above the base,
    /// [s_dn, s1] below with the same fractions; massive F is the identity outside [s0, s1].
    #[test]
    fn ends_apron_and_massive_periods() {
        let ex = |apron: (f64, f64), massive: bool| Some(Exact { u: [0.0, 0.6], f: [0.0, 0.45], apron, massive });
        let a = Corner { r: 0.75, q: 1.375, ladder: vec![20.0, 100.0, 180.0], t3: None, ex: ex((50.0, 100.0), false) };
        assert_eq!(a.period(128.0, 60.0), (50.0, 100.0));
        assert_eq!(a.period(128.0, 50.0), (50.0, 100.0));
        assert_eq!(a.period(128.0, 49.0), (20.0, 100.0));
        assert_eq!(a.period(128.0, 120.0), (100.0, 180.0));
        assert!((a.f(128.0, 50.0).0 - 50.0).abs() < 1e-9);
        assert!((a.f(128.0, 49.0).0 - (20.0 + 80.0 * 0.75 * 29.0 / 80.0)).abs() < 1e-9);
        assert_eq!(a.f(128.0, 99.0).1, 2);
        assert!(!a.has_tread());
        let m = Corner { r: 0.75, q: 1.375, ladder: vec![10.0, 90.0], t3: None, ex: ex((f64::NAN, f64::NAN), true) };
        assert_eq!(m.f(128.0, 5.0), (5.0, 1));
        assert_eq!(m.f(128.0, 95.0), (95.0, 2));
        assert!((m.f(128.0, 50.0).0 - (10.0 + 80.0 * 0.75 * 0.5)).abs() < 1e-9);
        assert_eq!(m.f_slope(128.0, 95.0), 1.0);
        assert!((m.f_slope(128.0, 50.0) - 0.75).abs() < 1e-9);
        let (p, r, q) = m.limb_slopes();
        assert!(p.is_nan() && (r - 0.75).abs() < 1e-9 && (q - 1.375).abs() < 1e-9);
    }

    /// Lane PERIOD-ENDS-2 (ct_sharp_ends_label): a hidden riser between two corners whose ladders disagree beyond twice the merge
    /// tolerance: where the lighter corner is on its cliff limb and the heavier on its talus, the heaviest-corner label says talus
    /// while the blended F climbs at the cliff-like slope; the label sub-key calls it cliff. Where both are on talus both say talus.
    #[test]
    fn ends_label_follows_the_blended_slope() {
        let n = 4;
        let ex = Some(Exact { u: [0.0, 0.6], f: [0.0, 0.3], apron: (f64::NAN, f64::NAN), massive: false });
        let lad = |o: f64| (0..6).map(|k| o + 100.0 * k as f64).collect::<Vec<f64>>();
        let corners = vec![Corner { r: 0.5, q: 1.75, ladder: lad(0.0), t3: None, ex }, Corner { r: 0.5, q: 1.75, ladder: lad(40.0), t3: None, ex }];
        let corner: Vec<u32> = (0..n * n).map(|i| if i % n < 2 { 0 } else { 1 }).collect();
        let mk = |lab: bool| SharpField { same_period: false, n, gx0: 0, gy0: 0, step_m: 100.0, d: vec![0.0; n * n], corner: corner.clone(), corners: corners.clone(), corner_ta: vec![1.5, 1.5], tread: false, stats: serde_json::Value::Null, blend: 1, tol: 6.0, soil_edge: 0.0, soil_tw: 0.0, soil: Vec::new(), ends_label: lab, dh_post: Vec::new(), dh_mode: 0 };
        let (a, b) = (mk(false), mk(true));
        // x 32 + 0.4 x 32: corner 0 (weight 0.6) at h2 190 is on its talus ([100, 200]: u 0.9 -> cliff!), pick h2 so that corner 0
        // (ladder 0, 100, ..) is on talus (u < 0.6) and corner 1 (40, 140, ..) on its cliff (u >= 0.6)
        let (x, y) = (32.0 + 0.4 * 32.0, 48.0);
        let h2 = 105.0; // corner 0: [100, 200] u 0.05 talus; corner 1: [40, 140] u 0.65 cliff
        let (_, ra) = a.apply(x, y, h2).unwrap();
        let (_, rb) = b.apply(x, y, h2).unwrap();
        assert_eq!(ra, 1);
        // blended slope 0.6 x 0.5 + 0.4 x 1.75 = 1.0 < sqrt(0.875) = 0.935? no: 1.0 > 0.935 -> cliff
        assert_eq!(rb, 2);
        let h2 = 120.0; // corner 1 u 0.8 cliff, corner 0 u 0.2 talus
        assert_eq!(b.apply(32.0 + 0.1 * 32.0, y, h2).unwrap().1, 1); // weight 0.9 on talus: 0.9 x 0.5 + 0.1 x 1.75 = 0.625
        let h2 = 110.0; // both on talus: corner 1 [40, 140] u 0.7 -> cliff; use 90: corner 0 [0,100] u 0.9 cliff; use 150
        let _ = h2;
        let (_, r1) = b.apply(x, y, 150.0).unwrap(); // corner 0 [100,200] u 0.5 talus; corner 1 [140, 240] u 0.1 talus
        assert_eq!(r1, 1);
        assert_eq!(limb_of_slope(0.2, (0.1, 0.5, 2.0), 1), 3);
        assert_eq!(limb_of_slope(0.5, (0.1, 0.5, 2.0), 3), 1);
        assert_eq!(limb_of_slope(1.5, (f64::NAN, 0.5, 2.0), 1), 2);
        assert_eq!(limb_of_slope(1.5, (f64::NAN, f64::NAN, 2.0), 1), 1);
    }

    /// Lane SPLIT-2M: ct_sharp_blend_period needs ct_sharp_blend=phase (given first) and is not serialized when off.
    #[test]
    fn blend_period_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_blend_period", "1").is_err(), "needs ct_sharp");
        q.set("ct_sharp", "1").unwrap();
        assert!(q.set("ct_sharp_blend_period", "1").is_err(), "needs ct_sharp_blend=phase");
        q.set("ct_sharp_blend", "phase").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("blend_period").is_none());
        q.set("ct_sharp_blend_period", "1").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"]["blend_period"].is_string());
        q.set("ct_sharp_blend_period", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("blend_period").is_none());
    }

    /// Lane SPLIT-2M (SHARP_BLEND_PERIOD_RULE): the measured ladder jump of the vs=1 split window (preview-client south wall, posts
    /// (1696, 4736) and (1728, 4736) of W9_vs5_s1 w_0_0, published ends and knots copied): A's period [-1347.4, 322.4] holds no
    /// riser at A; B has the extra level -109.4 (its bed is a riser). Both posts are on their talus limbs and their own staircase
    /// shifts are -43.6 and -34.3 m. Without the key A borrows B's knot fractions (B's own 431 m period) on its 1670 m period and
    /// the 2 m surface between them drops far below both shifts (the bead band); with it the change stays between the two own
    /// changes, the posts are reproduced and the row is continuous.
    #[test]
    fn blend_period_stops_the_ladder_jump_band() {
        let n = 4;
        let exa = Some(Exact { u: [0.0, (305.7 + 1347.5) / (322.3 + 1347.5)], f: [0.0, (259.1 + 1347.5) / (322.3 + 1347.5)], apron: (f64::NAN, f64::NAN), massive: false });
        let exb = Some(Exact { u: [0.0, (305.6 + 109.4) / (322.1 + 109.4)], f: [0.0, (258.9 + 109.4) / (322.1 + 109.4)], apron: (f64::NAN, f64::NAN), massive: false });
        let corners = vec![Corner { r: 0.9, q: 3.0, ladder: vec![-1480.8, -1347.5, 322.3, 432.1], t3: None, ex: exa },
                           Corner { r: 0.9, q: 3.0, ladder: vec![-1347.7, -109.4, 322.1, 431.9], t3: None, ex: exb }];
        let corner: Vec<u32> = (0..n * n).map(|i| if i % n < 2 { 0 } else { 1 }).collect();
        let mk = |sp: bool| SharpField { same_period: sp, n, gx0: 0, gy0: 0, step_m: 128.0, d: vec![0.0; n * n], corner: corner.clone(), corners: corners.clone(), corner_ta: vec![1.5, 1.5], tread: false, stats: serde_json::Value::Null, blend: 1, tol: 6.0, soil_edge: 0.0, soil_tw: 0.0, soil: Vec::new(), ends_label: false, dh_post: Vec::new(), dh_mode: 0 };
        let (off, on) = (mk(false), mk(true));
        let (ha, hb) = (198.4, 196.2);
        let own = |k: usize, h: f64| corners[k].f(128.0, h).0 - h;
        let (da, db) = (own(0, ha), own(1, hb));
        assert!((da + 43.6).abs() < 0.3 && (db + 34.3).abs() < 0.3, "own shifts {da} {db}");
        // the T1b posts are reproduced with or without the key
        for f in [&off, &on] {
            assert!((f.apply(32.0, 48.0, ha).unwrap().0 - (ha + da)).abs() < 1e-9);
            assert!((f.apply(64.0, 48.0, hb).unwrap().0 - (hb + db)).abs() < 1e-9);
        }
        // mid-cell: without the key the change falls far below both own changes; with it it stays between them (+- the 2.2 m
        // height difference x the talus F-slope deviation)
        let h2 = 0.5 * (ha + hb);
        let d_off = off.apply(48.0, 48.0, h2).unwrap().0 - h2;
        let d_on = on.apply(48.0, 48.0, h2).unwrap().0 - h2;
        assert!(d_off < da.min(db) - 10.0, "key off shows the band: {d_off}");
        assert!(d_on > da.min(db) - 0.5 && d_on < da.max(db) + 0.5, "key on between the own changes: {d_on} ({da}, {db})");
        // continuous along the row (h2 linear between the posts), no step larger than the own-shift gap per 0.25 m
        let mut prev = on.apply(32.0, 48.0, ha).unwrap().0;
        for k in 1..=128 {
            let t = k as f64 / 128.0;
            let z = on.apply(32.0 + 32.0 * t, 48.0, ha + (hb - ha) * t).unwrap().0;
            assert!((z - prev).abs() < 0.5, "step {} at t {t}", z - prev);
            prev = z;
        }
    }

    /// Lane CELL-RISER: ct_sharp_cell_riser needs ct_sharp_blend=phase and ct_sharp_ends=1 (given first), is not serialized when
    /// off, its D sub-key needs it, and ct_sharp_ends=0 drops both.
    #[test]
    fn cell_riser_key_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        assert!(q.set("ct_sharp_cell_riser", "1").is_err());
        q.set("ct_sharp", "1").unwrap();
        assert!(q.set("ct_sharp_cell_riser", "1").is_err());
        q.set("ct_sharp_ends", "1").unwrap();
        assert!(q.set("ct_sharp_cell_riser", "1").is_err());
        q.set("ct_sharp_blend", "phase").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("cell_riser").is_none());
        assert!(q.set("ct_sharp_cell_riser_m", "10").is_err());
        q.set("ct_sharp_cell_riser", "0").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("cell_riser").is_none());
        q.set("ct_sharp_cell_riser", "1").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["cell_riser"], serde_json::json!(SHARP_CELL_RISER_RULE));
        assert!(q.set("ct_sharp_cell_riser_m", "0").is_err());
        q.set("ct_sharp_cell_riser_m", "12.5").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["cell_riser_m"], serde_json::json!(12.5));
        assert!(q.set("ct_sharp_cell_riser_dh", "nearest").is_err());
        q.set("ct_sharp_cell_riser_dh", "smooth").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("cell_riser_dh").is_none());
        q.set("ct_sharp_cell_riser_dh", "bicubic").unwrap();
        assert_eq!(serde_json::to_value(&q).unwrap()["ct_sharp"]["cell_riser_dh"], serde_json::json!("bicubic"));
        assert_eq!(cell_riser_dh_mode(None), 2);
        q.set("ct_sharp_ends", "0").unwrap();
        let v = serde_json::to_value(&q).unwrap();
        assert!(v["ct_sharp"].get("cell_riser").is_none() && v["ct_sharp"].get("cell_riser_m").is_none() && v["ct_sharp"].get("cell_riser_dh").is_none());
        assert_eq!(SharpParams::default().cell_riser, None);
    }

    /// Lane CELL-RISER test field: h = 10 + g x (one period [0, 200] m, ladder -200 .. 400), every post an exact 2-limb corner
    /// whose knots depend on the post: `pick(i, j)` true = the "buried" profile (cliff foot f2 0.8 at u2 0.9, talus F-slope 0.89)
    /// else the "exposed" one (f2 0.3 at u2 0.5, F-slope 0.6): neighbouring posts at nearly the same base height then sit tens of
    /// metres apart on their own talus limbs (the s2_09 case). Returns h, z, regime and the ends.
    fn riser_field(n: usize, g: f64, pick: &dyn Fn(usize, usize) -> bool) -> (Vec<f64>, Vec<f64>, Vec<u8>, crate::input::CtEnds) {
        let nn = n * n;
        let (mut h, mut z, mut rg) = (vec![0.0; nn], vec![0.0; nn], vec![0u8; nn]);
        let mm = |v: f64| (v * 1000.0).round() as i32;
        let mut e = crate::input::CtEnds { s0: vec![0; nn], s1: vec![mm(200.0); nn], s_dn: vec![mm(-200.0); nn], s_up: vec![mm(400.0); nn], h1: vec![0; nn], z1: vec![0; nn], h2: vec![0; nn], z2: vec![0; nn], grad_permille: vec![(g * 1000.0).round() as u16; nn], kind: vec![1; nn] };
        for j in 0..n { for i in 0..n {
            let k = j * n + i;
            let (u2, f2) = if pick(i, j) { (0.9, 0.8) } else { (0.5, 0.3) };
            let hh = 10.0 + g * i as f64 * crate::T1B_M;
            let u = hh / 200.0;
            let (zz, r) = if u < u2 { (200.0 * u * f2 / u2, 1) } else { (200.0 * (f2 + (u - u2) * (1.0 - f2) / (1.0 - u2)), 2) };
            h[k] = (hh * 1000.0).round() / 1000.0; z[k] = (zz * 1000.0).round() / 1000.0; rg[k] = r;
            e.h2[k] = mm(200.0 * u2); e.z2[k] = mm(200.0 * f2);
        } }
        (h, z, rg, e)
    }

    /// Lane CELL-RISER: on the s2_09-like field (buried / exposed profiles alternating along the strike-parallel y axis) the key
    /// re-profiles the buried posts on the exposed neighbour's knots (dh smoothstep, the default; bilinear / bicubic printed): T1b
    /// posts are reproduced (< 5 mm), the ramps between the rows (talus-labelled 2 m slopes far steeper than the talus F-slope allows) become talus + one cliff limb, and the surface is
    /// continuous (max step halves when the sampling step halves: no jumps at cell edges or mid-lines).
    #[test]
    fn cell_riser_turns_the_ramp_into_talus_and_one_riser() {
        let n = 10;
        let g = 0.65;
        let (h, z, rg, e) = riser_field(n, g, &|_, j| j % 2 == 1);
        let c = CtParamsV1::v1();
        let (cp, lit) = (vec![500u16; n * n], vec![3u8; n * n]);
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.ends = Some(&e);
        let mut pe = phase_params();
        pe.ends = Some(SHARP_ENDS_RULE.to_string());
        let off = SharpField::build(&ci, &pe);
        let mut pc = pe.clone();
        pc.cell_riser = Some(SHARP_CELL_RISER_RULE.to_string());
        let on = SharpField::build(&ci, &pc);
        assert_eq!(on.dh_mode, 2);
        assert!(off.stats.get("cell_riser").is_none() && off.dh_post.is_empty());
        let st = &on.stats["cell_riser"];
        eprintln!("cell riser stats {st}");
        assert_eq!(st["reprofiled_posts"], serde_json::json!(n * n / 2));
        // T1b posts reproduced
        let z2at = |sf: &SharpField, x: f64, y: f64| 10.0 + g * x + sf.shift_bicubic(x, y);
        for j in 2..n - 2 { for i in 2..n - 2 {
            let k = j * n + i;
            let (x, y) = (i as f64 * crate::T1B_M, j as f64 * crate::T1B_M);
            let (zn, _) = on.apply(x, y, z2at(&on, x, y)).unwrap();
            assert!((zn - z[k]).abs() < 0.005, "post {i},{j}: {zn} vs {}", z[k]);
            let (zo, _) = off.apply(x, y, z2at(&off, x, y)).unwrap();
            assert!((zo - z[k]).abs() < 0.005, "key off post {i},{j}: {zo} vs {}", z[k]);
        } }
        // 2 m slopes along y (the strike: h is constant along y) over the interior: talus-labelled steepest, cliff share
        let prof = |sf: &SharpField, st: f64| -> (f64, f64, usize, f64) {
            let (mut tal_max, mut steps_max, mut cl, mut tot) = (0.0f64, 0.0f64, 0usize, 0usize);
            let x0 = 2.0 * crate::T1B_M;
            let mut xs = x0;
            while xs < (n - 3) as f64 * crate::T1B_M {
                let mut y = 2.0 * crate::T1B_M;
                let mut prev: Option<f64> = None;
                while y < (n - 3) as f64 * crate::T1B_M {
                    let (zz, r) = sf.apply(xs, y, z2at(sf, xs, y)).unwrap();
                    let (z1, _) = sf.apply(xs, y + st, z2at(sf, xs, y + st)).unwrap();
                    let s = ((z1 - zz) / st).abs();
                    if r == 1 { tal_max = tal_max.max(s); }
                    if r == 2 { cl += 1; }
                    tot += 1;
                    if let Some(p) = prev { steps_max = steps_max.max((zz - p).abs()); }
                    prev = Some(zz);
                    y += st;
                }
                xs += 7.0;
            }
            (tal_max, steps_max, cl, cl as f64 / tot as f64)
        };
        let (t_off, _, _, c_off) = prof(&off, 2.0);
        let (t_on, s2, _, c_on) = prof(&on, 2.0);
        let (_, s1, _, _) = prof(&on, 1.0);
        let (_, s05, _, _) = prof(&on, 0.5);
        eprintln!("talus-labelled max 2 m slope along strike: off {:.1} deg, on {:.1} deg; cliff share off {:.3} on {:.3}; max step 2/1/0.5 m: {s2:.2} {s1:.2} {s05:.2}",
            t_off.atan().to_degrees(), t_on.atan().to_degrees(), c_off, c_on);
        let mut t_bil = 0.0;
        for (m, nm) in [(0u8, "bilinear"), (1, "bicubic")] {
            let o2 = SharpField::build(&ci, &SharpParams { cell_riser_dh: Some(nm.to_string()), ..pc.clone() });
            assert_eq!(o2.dh_mode, m);
            let (t, s2m, _, c) = prof(&o2, 2.0);
            let (_, s05m, _, _) = prof(&o2, 0.5);
            if m == 0 { t_bil = t; }
            eprintln!("  dh {nm}: talus-labelled max {:.1} deg, cliff share {c:.3}, max step 2/0.5 m {s2m:.2} {s05m:.2}", t.atan().to_degrees());
        }
        // alternating shifts are the smoothstep's worst case (dh gradient 1.5x the bilinear one mid-cell): it still lowers the
        // steepest talus; the bilinear blend brings it to talus angles
        assert!(t_off.atan().to_degrees() > 55.0, "the ramp exists without the key");
        assert!(t_on < t_off, "talus-labelled slopes with the key (smooth): {}", t_on.atan().to_degrees());
        assert!(t_bil.atan().to_degrees() < 45.0, "talus-labelled slopes with the key (bilinear): {}", t_bil.atan().to_degrees());
        assert!(c_on > c_off, "a cliff limb appears");
        assert!(s05 < 0.6 * s1 && s1 < 0.6 * s2 || s05 < 0.35, "continuity: {s2} {s1} {s05}");
        // monotone F per post and the inverse
        for ck in &on.corners {
            let mut last = f64::NEG_INFINITY;
            for k in 0..400 { let hh = -150.0 + k as f64; let (v, _) = ck.f(on.step_m, hh); assert!(v >= last - 1e-9); last = v;
                assert!((ck.f_inv(on.step_m, v) - hh).abs() < 1e-6 || ck.f_slope(on.step_m, hh) == 0.0); }
        }
    }

    /// Lane CELL-RISER: posts whose foot is within D of the lowest same-period neighbour, the lowest post itself, apron / massive
    /// posts and posts whose neighbours sit on another period are untouched (own knots, dh 0); with no re-profiled post the field
    /// is the key-off field at every 2 m post (bit for bit).
    #[test]
    fn cell_riser_leaves_agreeing_and_other_period_posts() {
        let n = 8;
        let g = 0.65;
        let (h, z, rg, mut e) = riser_field(n, g, &|i, _| i >= 4);
        let c = CtParamsV1::v1();
        let (cp, lit) = (vec![500u16; n * n], vec![3u8; n * n]);
        // columns 4.. on another period (s0 shifted by 30 m): their neighbours in column 3 do not count
        for j in 0..n { for i in 4..n { e.s0[j * n + i] += 30_000; } }
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.ends = Some(&e);
        let mut pe = phase_params();
        pe.ends = Some(SHARP_ENDS_RULE.to_string());
        let off = SharpField::build(&ci, &pe);
        let mut pc = pe.clone();
        pc.cell_riser = Some(SHARP_CELL_RISER_RULE.to_string());
        let on = SharpField::build(&ci, &pc);
        assert_eq!(on.stats["cell_riser"]["reprofiled_posts"], serde_json::json!(0));
        assert!(on.dh_post.iter().all(|v| *v == 0.0));
        for k in 0..60 { for l in 0..60 {
            let (x, y) = (k as f64 * 3.7, l as f64 * 3.3);
            let z2 = 10.0 + g * x + off.shift_bicubic(x, y);
            assert_eq!(on.apply_ta(x, y, z2), off.apply_ta(x, y, z2));
        } }
    }
    /// Lane CELL-RISER-2: the donor / cap / monotone sub-keys need ct_sharp_cell_riser=1 (given first), are not serialized when
    /// at their defaults, and ct_sharp_cell_riser=0 drops them.
    #[test]
    fn cell_riser2_keys_opt_in() {
        let mut q = crate::DetailParamsV1::v3();
        q.set("ct_sharp", "1").unwrap();
        q.set("ct_sharp_blend", "phase").unwrap();
        q.set("ct_sharp_ends", "1").unwrap();
        assert!(q.set("ct_sharp_cell_riser_donor", "median").is_err());
        assert!(q.set("ct_sharp_cell_riser_cap", "1").is_err());
        q.set("ct_sharp_cell_riser", "1").unwrap();
        assert!(q.set("ct_sharp_cell_riser_donor", "lowest").is_err());
        assert!(q.set("ct_sharp_cell_riser_cap", "0").is_err());
        assert!(q.set("ct_sharp_cell_riser_cap", "9").is_err());
        q.set("ct_sharp_cell_riser_donor", "min").unwrap();
        assert!(serde_json::to_value(&q).unwrap()["ct_sharp"].get("cell_riser_donor").is_none());
        q.set("ct_sharp_cell_riser_donor", "median").unwrap();
        q.set("ct_sharp_cell_riser_cap", "0.5").unwrap();
        q.set("ct_sharp_cell_riser_dh", "monotone").unwrap();
        let v = serde_json::to_value(&q).unwrap();
        assert_eq!(v["ct_sharp"]["cell_riser_donor"], serde_json::json!("median"));
        assert_eq!(v["ct_sharp"]["cell_riser_cap"], serde_json::json!(0.5));
        assert_eq!(v["ct_sharp"]["cell_riser_dh"], serde_json::json!("monotone"));
        assert_eq!(cell_riser_dh_mode(Some("monotone")), 3);
        q.set("ct_sharp_cell_riser", "0").unwrap();
        let v = serde_json::to_value(&q).unwrap();
        assert!(v["ct_sharp"].get("cell_riser_donor").is_none() && v["ct_sharp"].get("cell_riser_cap").is_none());
        let d = SharpParams::default();
        assert!(d.cell_riser_donor.is_none() && d.cell_riser_cap.is_none());
    }

    /// Lane CELL-RISER-2: the monotone cubic passes through the posts, never leaves the two posts' range, is monotone where the
    /// data are, and has a continuous slope across the knot (one-sided difference quotients agree).
    #[test]
    fn cell_riser2_mono_cubic() {
        let cases = [[0.0, 0.0, 10.0, 10.0], [0.0, 10.0, 0.0, 10.0], [-5.0, 0.0, 40.0, 41.0], [3.0, 3.0, 3.0, 3.0], [0.0, -2.0, -60.0, -61.0]];
        for v in cases {
            assert!((mono_cubic(v, 0.0) - v[1]).abs() < 1e-12 && (mono_cubic(v, 1.0) - v[2]).abs() < 1e-12);
            let (lo, hi) = (v[1].min(v[2]), v[1].max(v[2]));
            let mut last = mono_cubic(v, 0.0);
            for k in 1..=100 {
                let a = mono_cubic(v, k as f64 / 100.0);
                assert!(a >= lo - 1e-9 && a <= hi + 1e-9, "{v:?} overshoot {a}");
                assert!((a - last) * (v[2] - v[1]) >= -1e-12, "{v:?} not monotone");
                last = a;
            }
        }
        // C1 at a knot: the slope leaving [v0..v3] at t = 1 equals the slope entering the next interval at t = 0
        let d = [0.0, 1.0, 5.0, 6.0, 20.0];
        let eps = 1e-6;
        let left = (mono_cubic([d[0], d[1], d[2], d[3]], 1.0) - mono_cubic([d[0], d[1], d[2], d[3]], 1.0 - eps)) / eps;
        let right = (mono_cubic([d[1], d[2], d[3], d[4]], eps) - mono_cubic([d[1], d[2], d[3], d[4]], 0.0)) / eps;
        assert!((left - right).abs() < 1e-4, "{left} {right}");
    }

    /// Lane CELL-RISER-2: on the alternating s2_09-like field, median / pair donors re-profile the same buried rows as the lowest
    /// foot (every neighbour row agrees); a single low outlier is the lowest foot's donor but not the median's / pair's; the cap
    /// drops posts whose shift would move the break by more than cap cells; dh=monotone reproduces the T1b posts and stays
    /// within the cell's corner shifts (no overshoot).
    #[test]
    fn cell_riser2_donor_cap_monotone() {
        let n = 10;
        let g = 0.65;
        let (h, z, rg, e) = riser_field(n, g, &|_, j| j % 2 == 1);
        let c = CtParamsV1::v1();
        let (cp, lit) = (vec![500u16; n * n], vec![3u8; n * n]);
        let mut ci = inputs(n, &h, &z, &rg, &cp, &lit, &c);
        ci.ends = Some(&e);
        let mut pc = phase_params();
        pc.ends = Some(SHARP_ENDS_RULE.to_string());
        pc.cell_riser = Some(SHARP_CELL_RISER_RULE.to_string());
        let base = SharpField::build(&ci, &pc);
        assert!(base.stats["cell_riser"].get("v2").is_none());
        let n_base = base.stats["cell_riser"]["reprofiled_posts"].as_u64().unwrap();
        // pair: the alternating rows confirm each other (two exposed neighbours in every 8-neighbourhood); median: a buried post has
        // 3 exposed + 2 buried neighbours in the interior (column edges fewer), the lower median is exposed where >= half are
        for (nm, lo) in [("pair", n_base), ("median", 1)] {
            let f = SharpField::build(&ci, &SharpParams { cell_riser_donor: Some(nm.to_string()), ..pc.clone() });
            let st = &f.stats["cell_riser"];
            eprintln!("donor {nm}: {st}");
            assert_eq!(st["v2"]["donor"], serde_json::json!(nm));
            let k = st["reprofiled_posts"].as_u64().unwrap();
            assert!(k >= lo && k <= n_base, "{nm}: {k} of {n_base}");
        }
        // one low outlier: all posts on one profile except (4, 4) exposed (foot 0.3 x 200 = 60 m vs 160 m)
        let (h1, z1, rg1, e1) = riser_field(n, g, &|i, j| !(i == 4 && j == 4));
        let mut c1 = inputs(n, &h1, &z1, &rg1, &cp, &lit, &c);
        c1.ends = Some(&e1);
        let m1 = SharpField::build(&c1, &pc);
        assert_eq!(m1.stats["cell_riser"]["reprofiled_posts"], serde_json::json!(8), "the lowest foot re-profiles the outlier's 8 neighbours");
        for nm in ["median", "pair"] {
            let f = SharpField::build(&c1, &SharpParams { cell_riser_donor: Some(nm.to_string()), ..pc.clone() });
            assert_eq!(f.stats["cell_riser"]["reprofiled_posts"], serde_json::json!(0), "{nm} ignores a single low outlier");
            assert!(f.dh_post.iter().all(|v| *v == 0.0));
        }
        // cap: |dh| <= cap x 32 x S
        let mx = base.dh_post.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        let cap_all = mx / (crate::T1B_M * g) * 1.01;
        let fa = SharpField::build(&ci, &SharpParams { cell_riser_cap: Some(cap_all), ..pc.clone() });
        assert_eq!(fa.stats["cell_riser"]["reprofiled_posts"].as_u64().unwrap(), n_base);
        assert_eq!(fa.dh_post, base.dh_post);
        let fz = SharpField::build(&ci, &SharpParams { cell_riser_cap: Some(0.5 * cap_all), ..pc.clone() });
        let k = fz.stats["cell_riser"]["reprofiled_posts"].as_u64().unwrap();
        assert_eq!(k + fz.stats["cell_riser"]["v2"]["capped_posts"].as_u64().unwrap(), n_base);
        assert!(fz.dh_post.iter().all(|v| v.abs() <= 0.5 * cap_all * crate::T1B_M * g + 1e-9));
        // monotone dh: posts reproduced, dh within the corner shifts
        let fm = SharpField::build(&ci, &SharpParams { cell_riser_dh: Some("monotone".into()), ..pc.clone() });
        assert_eq!(fm.dh_mode, 3);
        let z2at = |sf: &SharpField, x: f64, y: f64| 10.0 + g * x + sf.shift_bicubic(x, y);
        for j in 2..n - 2 { for i in 2..n - 2 {
            let (x, y) = (i as f64 * crate::T1B_M, j as f64 * crate::T1B_M);
            let (zn, _) = fm.apply(x, y, z2at(&fm, x, y)).unwrap();
            assert!((zn - z[j * n + i]).abs() < 0.005, "post {i},{j}: {zn} vs {}", z[j * n + i]);
        } }
        for k in 0..200 {
            let (x, y) = (70.0 + 1.03 * k as f64, 66.0 + 0.97 * k as f64);
            let (ix, iy) = ((x / crate::T1B_M).floor() as usize, (y / crate::T1B_M).floor() as usize);
            let cs = [fm.dh_post[iy * n + ix], fm.dh_post[iy * n + ix + 1], fm.dh_post[(iy + 1) * n + ix], fm.dh_post[(iy + 1) * n + ix + 1]];
            let (lo, hi) = cs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, b| (a.0.min(*b), a.1.max(*b)));
            let v = fm.dh_at(x, y);
            assert!(v >= lo - 1e-9 && v <= hi + 1e-9, "dh {v} outside [{lo}, {hi}]");
        }
    }
}
