//! Read-only VS3 observer over native carrier mutations. Existing codecs stay unchanged.
//! Temporal batches use the native homogeneous mixture rule within each layer.
//! This is an explicit conditional tracer extension, not recovered eruption dates.
use crate::layer_evolution::LayerChange;
use crate::layer_provenance::{Recorder, StackProvenance};
use crate::stratigraphy::{ParcelAmounts, StrataLayer, StrataParams, StrataStack, StrataV2Params};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
type Result<T> = std::result::Result<T, String>;
fn number(v: &Value) -> Result<f64> {
    let s = v.as_str().ok_or("history float must be encoded binary64")?;
    let x = f64::from_bits(u64::from_str_radix(s, 16).map_err(|_| "history float encoding")?);
    if !x.is_finite() {
        return Err("nonfinite history value".into());
    }
    Ok(x)
}
fn unsigned(v: &Value) -> Result<u64> {
    v.as_u64().ok_or("expected native unsigned field".into())
}
fn boolean(v: &Value) -> Result<bool> {
    v.as_bool().ok_or("expected native boolean field".into())
}
fn stack(v: &Value) -> Result<StrataStack> {
    Ok(StrataStack {
        basement_m: number(&v["basement_m"])?,
        basement_lith: unsigned(&v["basement_lith"])? as u8,
        basement_age_my: number(&v["basement_age_my"])?,
        last_exposure_my: number(&v["last_exposure_my"])?,
        rock_uplift_m: number(&v["rock_uplift_m"])?,
        elev_ref_m: number(&v["elev_ref_m"])?,
        cover: v["cover"]
            .as_array()
            .ok_or("native cover array")?
            .iter()
            .map(|l| {
                Ok(StrataLayer {
                    thickness_m: number(&l["thickness_m"])?,
                    age_base_my: number(&l["age_base_my"])?,
                    age_top_my: number(&l["age_top_my"])?,
                    max_burial_m: number(&l["max_burial_m"])?,
                    lith: unsigned(&l["lith"])? as u8,
                    env: unsigned(&l["env"])? as u8,
                    grade: unsigned(&l["grade"])? as u8,
                    flags: unsigned(&l["flags"])? as u8,
                })
            })
            .collect::<Result<Vec<_>>>()?,
    })
}
fn amounts(v: &Value) -> Result<ParcelAmounts> {
    Ok(ParcelAmounts {
        erosion: number(&v["erosion"])?,
        deposit: number(&v["deposit"])?,
        class: unsigned(&v["class"])? as u8,
        arc: number(&v["arc"])?,
        plume: number(&v["plume"])?,
        thinning: number(&v["thinning"])?,
        basement_add: number(&v["basement_add"])?,
        short_pos: number(&v["short_pos"])?,
        short_neg: number(&v["short_neg"])?,
        belt_rest: number(&v["belt_rest"])?,
        land: boolean(&v["land"])?,
        warm_shelf: boolean(&v["warm_shelf"])?,
    })
}
// Native reconciliation may temporarily carry a negative basement deficit.
// The unchanged provenance observer represents that as zero material. Preserve
// the signed native stack; validate only its material quantity against max(0,b).
fn validate_observed_material(p: &StackProvenance, st: &StrataStack) -> Result<()> {
    p.basement.validate(st.basement_m.max(0.))?;
    if p.cover.len() != st.cover.len() {
        return Err("observed cover shape mismatch".into());
    }
    for (a, b) in p.cover.iter().zip(&st.cover) {
        a.validate(b.thickness_m.max(0.))?;
    }
    Ok(())
}

fn params(v: &Value) -> Result<StrataParams> {
    macro_rules! fields { ($ty:ident, $v:expr, $($f:ident),+) => {$ty {$($f:number(&$v[stringify!($f)])?,)+}}; }
    let w = &v["v2"];
    Ok(StrataParams {
        max_layers: unsigned(&v["max_layers"])? as usize,
        v2: if w.is_null() {
            None
        } else {
            Some(fields!(
                StrataV2Params,
                w,
                shelf_width_m,
                slope_foot_depth_m,
                shelf_majority_fraction,
                platform_cover_max_m,
                platform_full_elev_m,
                platform_zero_elev_m,
                platform_sandstone_fraction,
                platform_carbonate_fraction,
                platform_age_base_my,
                platform_age_top_my
            ))
        },
        merge_gap_my: number(&v["merge_gap_my"])?,
        unconformity_gap_my: number(&v["unconformity_gap_my"])?,
        min_layer_m: number(&v["min_layer_m"])?,
        trace_m: number(&v["trace_m"])?,
        arc_extrusive_fraction: number(&v["arc_extrusive_fraction"])?,
        plume_extrusive_fraction: number(&v["plume_extrusive_fraction"])?,
        shortening_pluton_fraction: number(&v["shortening_pluton_fraction"])?,
        pluton_roof_m: number(&v["pluton_roof_m"])?,
        metamorphic_depth_m: number(&v["metamorphic_depth_m"])?,
        high_grade_depth_m: number(&v["high_grade_depth_m"])?,
        shelf_depth_m: number(&v["shelf_depth_m"])?,
        carbonate_max_paleolat_deg: number(&v["carbonate_max_paleolat_deg"])?,
        carbonate_max_supply_m_per_my: number(&v["carbonate_max_supply_m_per_my"])?,
        arid_lat_min_deg: number(&v["arid_lat_min_deg"])?,
        arid_lat_max_deg: number(&v["arid_lat_max_deg"])?,
        humid_lat_max_deg: number(&v["humid_lat_max_deg"])?,
        humid_midlat_min_deg: number(&v["humid_midlat_min_deg"])?,
        humid_midlat_max_deg: number(&v["humid_midlat_max_deg"])?,
        coal_max_elevation_m: number(&v["coal_max_elevation_m"])?,
        conglomerate_relief_m: number(&v["conglomerate_relief_m"])?,
        river_mouth_area_m2: number(&v["river_mouth_area_m2"])?,
        publish_depth_m: number(&v["publish_depth_m"])?,
        basement_min_m: number(&v["basement_min_m"])?,
    })
}
/// Quantization of an observed native dry-mass proxy, same unit/density as VS2.
/// Integer identities are exact; the float input is NOT accurate to one ug.
pub fn mass(volume: f64, rho: f64) -> Result<i128> {
    let q = volume * rho * 1e9;
    if !q.is_finite() || q < 0. || q >= i128::MAX as f64 {
        return Err("VS3 mass range".into());
    }
    Ok(q.round_ties_even() as i128)
}
fn heights(p: &StackProvenance) -> BTreeMap<String, f64> {
    let mut h = BTreeMap::new();
    for c in std::iter::once(&p.basement)
        .chain(&p.cover)
        .flat_map(|p| &p.contributions)
    {
        *h.entry(c.episode.clone()).or_default() += c.thickness_m;
    }
    h
}
#[derive(Default)]
struct Account {
    initial: i128,
    emplaced: i128,
    kinematic_in: i128,
    eroded: i128,
    moved: i128,
    numeric_in: i128,
    numeric_out: i128,
    remaining: i128,
    peak: f64,
    meta: Value,
}
fn q(n: i128) -> Value {
    json!({"unit":"ug","version":1,"decimal":n.to_string()})
}
/// The producer authenticates the archive before calling this view. No source
/// identity is changed and no legacy decoder is tricked into accepting it.
pub fn observe(archive: &Value, rho: f64) -> Result<Value> {
    let p = params(&archive["params"])?;
    p.validate()?;
    let events = archive["events"].as_array().ok_or("events")?;
    let selected = archive["config"]["selected"].as_array().ok_or("selected")?;
    let mut rows = Vec::new();
    for id in selected {
        let id = unsigned(id)?;
        let es: Vec<_> = events.iter().filter(|e| e["id"] == id).collect();
        let first = *es.first().ok_or("missing birth")?;
        if !matches!(unsigned(&first["kind"])?, 0 | 2) {
            return Err("missing initial/born carrier".into());
        }
        let birth = number(&first["end_my"])?;
        let area = number(&first["area_m2"])?;
        if area <= 0. {
            return Err("invalid carrier area".into());
        }
        let mut last = stack(&first["after"])?;
        let mut provenance = crate::layer_provenance::initial(&last, id);
        let mut accounts = BTreeMap::<String, Account>::new();
        for c in std::iter::once(&mut provenance.basement)
            .chain(&mut provenance.cover)
            .flat_map(|p| &mut p.contributions)
        {
            c.episode = format!("imported-at-birth:{id}:{}", c.episode);
        }
        for (k, h) in heights(&provenance) {
            let n = mass(h * area, rho)?;
            accounts.insert(k,Account {initial:n,remaining:n,meta:json!({"origin":"initial_or_imported_crust_not_new_eruption","birth_my":birth}),..Default::default()});
        }
        let mut time = birth;
        let mut operations = 0;
        let mut changes = 0;
        for e in es.into_iter().skip(1) {
            if e["kind"] != 1 {
                return Err("retired/reborn carrier not supported by surviving source book".into());
            }
            let t0 = number(&e["start_my"])?;
            let t1 = number(&e["end_my"])?;
            if t0 != time || area != number(&e["area_m2"])? || stack(&e["before"])? != last {
                return Err("source history continuity mismatch".into());
            }
            let mut rec = Recorder {
                state: provenance,
                carrier: id,
                amounts: amounts(&e["amounts"])?,
                params: p,
                t0,
                t1,
            };
            for (index, c) in e["evolution"]["changes"]
                .as_array()
                .ok_or("missing native mutation trace")?
                .iter()
                .enumerate()
            {
                let kind = u8::try_from(unsigned(&c["kind"])?).map_err(|_| "kind range")?;
                if kind as usize >= crate::layer_evolution::NAMES.len() {
                    return Err("unknown native mutation".into());
                }
                let change = LayerChange {
                    kind,
                    links: serde_json::from_value(c["links"].clone()).map_err(|e| e.to_string())?,
                    before: stack(&c["before"])?,
                    after: stack(&c["after"])?,
                };
                validate_observed_material(&rec.state, &change.before)?;
                let before = heights(&rec.state);
                rec.capture(&change);
                // Prevent the native interval compressor from coalescing these
                // explicitly conditional tracer batches across different steps.
                for x in std::iter::once(&mut rec.state.basement)
                    .chain(&mut rec.state.cover)
                    .flat_map(|l| &mut l.contributions)
                {
                    if !before.contains_key(&x.episode) {
                        let source = x.episode.clone();
                        x.episode = format!("batch:{id}:{:016x}:{index}:{source}", t0.to_bits());
                        let volcanic = source.contains(":arc-volcanic:")
                            || source.contains(":plume-volcanic:");
                        let branch = if source.contains(":arc-") {
                            Some("arc")
                        } else if source.contains(":plume-") {
                            Some("plume")
                        } else {
                            None
                        };
                        accounts.entry(x.episode.clone()).or_insert_with(||Account {meta:json!({"origin":"native_addition","native_episode":source,"interval_my":[t0,t1],"branch":branch,"extrusive":volcanic,"kind":crate::layer_evolution::NAMES[kind as usize],"composition":null,"chemistry_status":"not recorded in native source"}),..Default::default()});
                    }
                }
                validate_observed_material(&rec.state, &change.after)?;
                let after = heights(&rec.state);
                for key in before
                    .keys()
                    .chain(after.keys())
                    .cloned()
                    .collect::<BTreeSet<_>>()
                {
                    let bh = before.get(&key).copied().unwrap_or(0.);
                    let ah = after.get(&key).copied().unwrap_or(0.);
                    let b = mass(bh * area, rho)?;
                    let a = mass(ah * area, rho)?;
                    let ac = accounts.get_mut(&key).ok_or("unowned source")?;
                    if ac.remaining != b {
                        return Err("source mass continuity".into());
                    }
                    let delta = a - b;
                    // Transformations without a material source/sink can change
                    // floating summation order. Bound and retain this separately.
                    let numeric = matches!(kind, 3 | 7 | 16 | 17 | 18 | 19 | 20 | 21);
                    if numeric
                        && (ah - bh).abs() > 64. * f64::EPSILON * (ah.abs() + bh.abs()).max(1.)
                    {
                        return Err("nonconservative native rearrangement".into());
                    }
                    if delta >= 0 {
                        if numeric {
                            ac.numeric_in += delta;
                        } else if !before.contains_key(&key) {
                            ac.emplaced += delta;
                        } else {
                            ac.kinematic_in += delta;
                        }
                    } else if numeric {
                        ac.numeric_out -= delta;
                    } else if kind <= 2 {
                        ac.eroded -= delta;
                    } else {
                        ac.moved -= delta;
                    }
                    ac.remaining = a;
                    ac.peak = ac.peak.max(ah * area);
                }
                changes += 1;
            }
            provenance = rec.state;
            last = stack(&e["after"])?;
            time = t1;
            validate_observed_material(&provenance, &last)?;
            operations += 1;
        }
        let mut out = Vec::new();
        let mut volcanic = 0i128;
        let mut imported = 0i128;
        for (key, a) in accounts {
            if a.initial + a.emplaced + a.kinematic_in + a.numeric_in
                != a.eroded + a.moved + a.numeric_out + a.remaining
            {
                return Err("integer source closure".into());
            }
            if a.meta["extrusive"] == true {
                volcanic += a.remaining;
            }
            if a.initial > 0 {
                imported += a.remaining;
            }
            out.push(json!({"batch_id":key,"metadata":a.meta,"initial_imported":q(a.initial),"emplaced":q(a.emplaced),"kinematic_import":q(a.kinematic_in),"eroded":q(a.eroded),"redistributed_out":q(a.moved),"numerical_in":q(a.numeric_in),"numerical_out":q(a.numeric_out),"remaining":q(a.remaining),"closure_ug":"0","native_peak_m3":a.peak}));
        }
        rows.push(json!({"carrier":id,"birth_kind":first["kind"],"birth_my":birth,"end_my":time,"area_m2":area,"operations":operations,"mutations":changes,"batches":out,"remaining_volcanic":q(volcanic),"remaining_imported":q(imported),"remaining_volcanic_m3":volcanic as f64/(rho*1e9),"birth_inputs":first["birth_inputs"]}));
    }
    Ok(
        json!({"schema":"pf.born-carrier-source-book/1","carriers":rows,"density_kg_m3":rho,
        "retention_model":"conditional time-batch tracer using unchanged native homogeneous-within-layer mixture transport; no recovered eruption schedule or grain ancestry",
        "quantity_basis":"native area times thickness times pinned legacy crust density; source quantization nearest-even ug; exact integer closure does not imply microgram physical precision",
        "redistribution":"exports have unknown destination; no recipient is credited without a transport record",
        "authority":"observation only; not kernel stock admission or cone placement"}),
    )
}
