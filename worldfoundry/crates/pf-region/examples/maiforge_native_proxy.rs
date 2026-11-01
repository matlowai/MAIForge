//! Bounded render-only residual controls sampled directly from the native model.
//! This never changes terrain laws and does not bake a global fine grid.
use pf_hydro_catchment::{input::EndpointInput, NONE};
use pf_region::{
    planet_catchment::{load_published, EndpointHydrology, World, MODEL, RESOLUTION},
    planet_relief::WorldRelief,
};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
    time::Instant,
};
const MAX_TOTAL: usize = 80_000_000;
const MAX_CONTEXT: usize = 12_000_000;
fn hash(bytes: &[u8]) -> String {
    pf_hydro_catchment::input::sha256_bytes(bytes)
}
fn work(
    input: &EndpointInput,
    hydro: &EndpointHydrology,
    ny: usize,
    nx: usize,
) -> Result<Value, String> {
    let g = input.grid;
    let mut outlet = vec![NONE; g.len()];
    // Work planning only: identical receiving-body partition and two-cell halo.
    for c in 0..g.len() {
        let mut path = Vec::new();
        let mut i = c;
        while outlet[i] == NONE && hydro.body[i] == NONE && input.published_flow[i] != NONE {
            path.push(i);
            i = input.published_flow[i];
            if i >= g.len() || path.len() > g.len() {
                return Err("receiver cycle".into());
            }
        }
        let key = if outlet[i] == NONE { i } else { outlet[i] };
        outlet[i] = key;
        for c in path {
            outlet[c] = key
        }
    }
    let wet = |c: usize| {
        hydro.body[c] != NONE && input.elevation[c] <= hydro.records[&hydro.body[c]].water_surface_m
    };
    let rows: Vec<usize> = (0..ny)
        .map(|i| {
            (((i as f64 + 0.5) * g.ny as f64 / ny as f64 - 0.5).round() as usize).min(g.ny - 1)
        })
        .collect();
    let cols: Vec<usize> = (0..nx)
        .map(|i| {
            (((i as f64 + 0.5) * g.nx as f64 / nx as f64 - 0.5).round() as usize).min(g.nx - 1)
        })
        .collect();
    let mut contexts = BTreeMap::<usize, Vec<usize>>::new();
    let mut controls = Vec::new();
    for &r in &rows {
        for &c in &cols {
            let cell = r * g.nx + c;
            let index = controls.len();
            let key = if wet(cell) { None } else { Some(outlet[cell]) };
            if let Some(k) = key {
                contexts.entry(k).or_default().push(index)
            }
            controls.push(json!({"cell":cell,"row":r,"col":c,"latitudeDegrees":-90.+(r as f64+0.5)*180./g.ny as f64,"longitudeDegrees":-180.+(c as f64+0.5)*360./g.nx as f64,"latitudeSnapOffsetDegrees":(r as f64+0.5)*180./g.ny as f64-((index/nx) as f64+0.5)*180./ny as f64,"longitudeSnapOffsetDegrees":(c as f64+0.5)*360./g.nx as f64-((index%nx) as f64+0.5)*360./nx as f64,"outletContext":key,"nativeWetFloor":wet(cell)}));
        }
    }
    let mut poles = Vec::new();
    for r in [0, g.ny - 1] {
        let c = (r * g.nx..(r + 1) * g.nx)
            .find(|&c| wet(c))
            .unwrap_or(r * g.nx);
        let key = if wet(c) { None } else { Some(outlet[c]) };
        if let Some(k) = key {
            contexts.entry(k).or_default();
        }
        poles.push(json!({"latitudeDegrees":if r==0{-90.}else{90.},"longitudeDegrees":-180.+((c%g.nx) as f64+0.5)*360./g.nx as f64,"rowNodeCell":c,"outletContext":key}));
    }
    let mut plans = Vec::new();
    for (key, indices) in contexts {
        let core: BTreeSet<usize> = outlet
            .iter()
            .enumerate()
            .filter_map(|(i, &k)| (k == key).then_some(i))
            .collect();
        let mut included = core.clone();
        for _ in 0..2 {
            let before = included.clone();
            for c in before {
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        if let Some(j) = g.neighbor(c, dx, dy) {
                            included.insert(j);
                        }
                    }
                }
            }
        }
        plans.push(json!({"outletCell":key,"coreCells":core.len(),"haloCells":included.len()-core.len(),"physicalNodes":included.len()*RESOLUTION*RESOLUTION,"controlIndices":indices}));
    }
    let total: usize = plans
        .iter()
        .map(|v| v["physicalNodes"].as_u64().unwrap() as usize)
        .sum();
    let largest = plans
        .iter()
        .map(|v| v["physicalNodes"].as_u64().unwrap())
        .max()
        .unwrap_or(0);
    Ok(
        json!({"schema":"maiforge.native-proxy-preflight/v1","model":MODEL,"controlGrid":[ny,nx],"nativeGrid":[g.ny,g.nx],"controls":controls,"poles":poles,"contexts":plans,"contextCount":plans.len(),"plannedNativeNodes":total,"maximumContextNodes":largest,"withinBudget":total<=MAX_TOTAL&&largest<=MAX_CONTEXT as u64,"bounds":{"totalNativeNodes":MAX_TOTAL,"perContextNativeNodes":MAX_CONTEXT,"memoryGiB":10,"wallSeconds":300},"planning":"Unique canonical owner contexts, two native-cell halo, grouped before sampling. Submerged wet-floor controls require no process solve; exact coarse-cell centres avoid owner junctions."}),
    )
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 {
        return Err(
            "usage: maiforge_native_proxy WORLD_ROOT OUTPUT_DIR ROWS(16|8) plan|export".into(),
        );
    }
    let root = Path::new(&a[1]);
    let out = Path::new(&a[2]);
    let ny = a[3].parse::<usize>()?;
    if ![8, 16].contains(&ny) || !["plan", "export"].contains(&a[4].as_str()) {
        return Err("unsupported bounded request".into());
    }
    if a[4] == "export" && ny != 8 {
        return Err("only the bounded 8x16 export is supported".into());
    }
    if out.exists() {
        return Err("output exists; choose fresh directory".into());
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .canonicalize()?;
    if parent
        .join(out.file_name().ok_or("output name")?)
        .starts_with(root.canonicalize()?)
    {
        return Err("output inside input prohibited".into());
    }
    let start = Instant::now();
    let (source, input, hydro) = load_published(root).map_err(std::io::Error::other)?;
    let plan = work(&input, &hydro, ny, ny * 2).map_err(std::io::Error::other)?;
    std::fs::create_dir(out)?;
    std::fs::write(
        out.join("preflight.json"),
        serde_json::to_vec_pretty(&plan)?,
    )?;
    println!(
        "{}",
        json!({"controlGrid":plan["controlGrid"],"contexts":plan["contextCount"],"nodes":plan["plannedNativeNodes"],"largestContext":plan["maximumContextNodes"],"withinBudget":plan["withinBudget"],"preflightSeconds":start.elapsed().as_secs_f64()})
    );
    if a[4] == "plan" {
        return Ok(());
    }
    if plan["withinBudget"] != true {
        return Err("native proxy preflight exceeds bounded work; no controls sampled".into());
    }
    let world = World::from_endpoint(input, hydro, &source).map_err(std::io::Error::other)?;
    let datum = WorldRelief::new(source).map_err(std::io::Error::other)?;
    let provenance = world.provenance();
    let controls = plan["controls"].as_array().ok_or("controls")?;
    let mut values = vec![None; controls.len()];
    let mut max_replay = 0.0f64;
    let mut receipts = Vec::new();
    let mut order: Vec<usize> = controls
        .iter()
        .enumerate()
        .filter_map(|(i, c)| c["outletContext"].is_null().then_some(i))
        .collect();
    for context in plan["contexts"].as_array().ok_or("contexts")? {
        order.extend(
            context["controlIndices"]
                .as_array()
                .ok_or("indices")?
                .iter()
                .map(|i| i.as_u64().unwrap() as usize),
        );
    }
    for (ordinal, index) in order.into_iter().enumerate() {
        if start.elapsed().as_secs() > 290 {
            return Err("native proxy wall budget reached; no proxy manifest published".into());
        }
        let c = &controls[index];
        let lat = c["latitudeDegrees"].as_f64().ok_or("latitude")?;
        let lon = c["longitudeDegrees"].as_f64().ok_or("longitude")?;
        let native = world.sample(lat, lon).map_err(std::io::Error::other)?;
        let replay = world.sample(lat, lon).map_err(std::io::Error::other)?;
        max_replay = max_replay.max((native - replay).abs());
        let base = datum.datum(lat, lon);
        values[index] = Some((native, base, native - base));
        if !c["outletContext"].is_null() {
            let key = c["outletContext"].as_u64().unwrap();
            if !receipts.iter().any(|r: &Value| r["outletCell"] == key) {
                let cell = c["cell"].as_u64().unwrap() as u32;
                let r = world
                    .receipts(&BTreeSet::from([cell]))
                    .map_err(std::io::Error::other)?;
                let r = r
                    .as_array()
                    .and_then(|r| r.first())
                    .ok_or("missing native receipt")?;
                let planned = plan["contexts"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|v| v["outletCell"] == key)
                    .unwrap();
                if r["physical_nodes"] != planned["physicalNodes"] {
                    return Err("actual native context work differs from preflight".into());
                }
                receipts.push(json!({"outletCell":key,"physicalNodes":r["physical_nodes"],"coreCells":r["core_cells"],"haloCells":r["halo_cells"]}));
            }
        }
        eprintln!(
            "control {}/{} elapsed {:.3}s",
            ordinal + 1,
            controls.len(),
            start.elapsed().as_secs_f64()
        );
    }
    let mut pole_values = Vec::new();
    for p in plan["poles"].as_array().unwrap() {
        let lat = p["latitudeDegrees"].as_f64().unwrap();
        let lon = p["longitudeDegrees"].as_f64().unwrap();
        let native = world.sample(lat, lon).map_err(std::io::Error::other)?;
        let base = datum.datum(lat, lon);
        pole_values.push(json!({"latitudeDegrees":lat,"longitudeDegrees":lon,"nativeM":native,"datumM":base,"residualM":native-base}));
    }
    let mut bytes = Vec::new();
    let mut quant = 0.0f64;
    let mut truth = Vec::new();
    for (i, v) in values.into_iter().enumerate() {
        let (native, base, residual) = v.ok_or("incomplete controls")?;
        let q = residual as f32;
        if !q.is_finite() {
            return Err("nonfinite residual".into());
        }
        quant = quant.max((residual - f64::from(q)).abs());
        bytes.extend(q.to_le_bytes());
        truth.push(json!({"index":i,"nativeM":native,"datumM":base,"residualM":residual}));
    }
    let rows: Vec<usize> = controls
        .iter()
        .step_by(ny * 2)
        .map(|c| c["row"].as_u64().unwrap() as usize)
        .collect();
    let cols: Vec<usize> = controls[..ny * 2]
        .iter()
        .map(|c| c["col"].as_u64().unwrap() as usize)
        .collect();
    let mut manifest = json!({"schema":"maiforge.viewer-native-proxy/v1","model":MODEL,"source":provenance,"seed":datum.input.seed,"radiusM":datum.input.radius_m,"nativeGrid":[datum.input.rows,datum.input.cols],"controlGrid":[ny,ny*2],"rowIndices":rows,"colIndices":cols,"latitudeDegrees":controls.iter().step_by(ny*2).map(|c|c["latitudeDegrees"].clone()).collect::<Vec<_>>(),"longitudeDegrees":controls[..ny*2].iter().map(|c|c["longitudeDegrees"].clone()).collect::<Vec<_>>(),"grid":{"order":"south-first","layout":"row-major","sampling":"native coarse-cell centres nearest regular control targets; exact coordinates above"},"sourceResolutionM":datum.input.radius_m*std::f64::consts::PI/(datum.input.rows*RESOLUTION) as f64,"predictor":"native endpoint datum: cell-centred bilinear with native pole-row means","interpolation":"render-only: bilinear residual between the declared latitude/longitude controls, periodic longitude; outer control row blends linearly to declared pole residual","poleSamples":pole_values,"fields":{"residual":{"file":"residual.f32","type":"f32le","unit":"m","bytes":bytes.len(),"sha256":hash(&bytes),"maxQuantizationErrorM":quant}},"qualification":"Sparse render-only predictor of the same native solid model. Global interpolation error is unknown. Not simulation, collision, water-stage, material or finite-stock authority.","maximumCanonicalReplayDifferenceM":max_replay});
    manifest["id"] = format!("blake3:{}", blake3::hash(&serde_json::to_vec(&manifest)?)).into();
    std::fs::write(out.join("residual.f32"), bytes)?;
    std::fs::write(
        out.join("controls.json"),
        serde_json::to_vec_pretty(&truth)?,
    )?;
    std::fs::write(
        out.join("proxy.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    std::fs::write(
        out.join("execution.json"),
        serde_json::to_vec_pretty(
            &json!({"elapsedSeconds":start.elapsed().as_secs_f64(),"actualContexts":receipts,"maximumCanonicalReplayDifferenceM":max_replay}),
        )?,
    )?;
    Ok(())
}
