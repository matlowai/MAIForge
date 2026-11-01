//! Export a compact browser data source from verified native planet artifacts.
//! Only numeric fields and allowlisted metadata are published; original headers,
//! source paths, execution records and historical notes remain in the input.
use pf_contracts::artifact::{
    read_artifact_header, read_descriptor_bytes, verify_artifact, ArtifactOnDisk, Dtype,
};
use pf_contracts::domain::{HydrologyDomainV0, WaterBodyKindV0, WaterBodyV0};
use pf_contracts::terrain::StratigraphyDomainV0;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    path::{Component, Path, PathBuf},
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const MAX_CELLS: usize = 16_777_216;
const SCHEMA: &str = "maiforge.viewer-world/v1";

fn sha(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn finite(v: f64, what: &str) -> Result<f64> {
    if !v.is_finite() {
        return Err(format!("nonfinite {what}").into());
    }
    Ok(v)
}
fn number(v: &Value, key: &str) -> Result<f64> {
    finite(
        v[key]
            .as_f64()
            .ok_or_else(|| format!("missing numeric {key}"))?,
        key,
    )
}
fn safe_relative(path: &str) -> Result<()> {
    if path.is_empty()
        || !Path::new(path)
            .components()
            .all(|c| matches!(c, Component::Normal(_)))
    {
        return Err("artifact filename must be a normal relative path".into());
    }
    Ok(())
}
fn checked_artifact(path: &Path) -> Result<ArtifactOnDisk<Value>> {
    // Validate all file references before the native verifier reads blob bytes.
    let a = read_artifact_header::<Value>(path)?;
    for f in &a.header.semantic.fields {
        safe_relative(&f.file)?;
    }
    for f in &a.header.semantic.sidecars {
        safe_relative(&f.path)?;
    }
    let digest = verify_artifact::<Value>(path)?;
    if digest != a.digest {
        return Err("source artifact changed during verification".into());
    }
    Ok(a)
}
fn field_grid(
    a: &ArtifactOnDisk<Value>,
    name: &str,
    dtype: Dtype,
    unit: &str,
) -> Result<(usize, usize, f64)> {
    let f = a
        .header
        .semantic
        .fields
        .iter()
        .find(|f| f.name == name)
        .ok_or_else(|| format!("missing field {name}"))?;
    if f.dtype != dtype || f.units != unit || f.shape.len() != 2 || f.order != "C" {
        return Err(format!("unsupported field contract {name}").into());
    }
    let rows = usize::try_from(f.shape[0])?;
    let cols = usize::try_from(f.shape[1])?;
    let n = rows.checked_mul(cols).ok_or("field grid overflow")?;
    if rows < 2 || cols < 4 || n > MAX_CELLS {
        return Err("unsupported field dimensions or export memory bound".into());
    }
    let d = a
        .header
        .semantic
        .domains
        .iter()
        .find(|d| d.domain_id == f.domain)
        .ok_or("missing field domain")?;
    let coords = serde_json::to_value(&d.discretization.coordinate_system)?;
    if coords != "SphericalLatLon"
        || d.discretization.backend != "latlon-v0"
        || d.discretization.element_count != n as u64
        || d.datum.center != [0.0; 3]
    {
        return Err("expected a planet-centred native latitude/longitude domain".into());
    }
    let radius = finite(d.datum.radius_m, "radius")?;
    if radius <= 0.0 {
        return Err("planet radius must be positive".into());
    }
    Ok((rows, cols, radius))
}
fn check_grid(got: (usize, usize, f64), want: (usize, usize, f64)) -> Result<()> {
    if got != want {
        return Err("source field grid or radius mismatch".into());
    }
    Ok(())
}
fn require_inputs(a: &ArtifactOnDisk<Value>, digests: &[&str]) -> Result<()> {
    if digests
        .iter()
        .any(|d| !a.header.semantic.inputs.iter().any(|v| v == d))
    {
        return Err("source artifact does not pin the required upstream data".into());
    }
    Ok(())
}
fn f32_bytes(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}
fn native_bytes(
    a: &ArtifactOnDisk<Value>,
    path: &Path,
    name: &str,
    dtype: Dtype,
) -> Result<Vec<u8>> {
    let field = a
        .header
        .semantic
        .fields
        .iter()
        .find(|f| f.name == name && f.dtype == dtype)
        .ok_or_else(|| format!("missing field or wrong encoding: {name}"))?;
    // Bind every consumed byte to the already verified header, even if a caller
    // changes source files while this read-only export is running.
    Ok(read_descriptor_bytes(path, &a.header.semantic, field)?)
}
fn checked_floats(
    a: &ArtifactOnDisk<Value>,
    path: &Path,
    name: &str,
    n: usize,
) -> Result<Vec<f32>> {
    let bytes = native_bytes(a, path, name, Dtype::F32Le)?;
    let values: Vec<f32> = bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    if values.len() != n || values.iter().any(|v| !v.is_finite()) {
        return Err(format!("invalid finite field {name}").into());
    }
    Ok(values)
}
fn water_surface(
    elevation: &[f32],
    ids: &[u32],
    bodies: &[WaterBodyV0],
    sea: f64,
    radius: f64,
) -> Result<(Vec<f32>, Vec<u8>, f64, f64)> {
    if elevation.len() != ids.len() {
        return Err("water/height shape mismatch".into());
    }
    let records: BTreeMap<_, _> = bodies.iter().map(|b| (b.body_id, b)).collect();
    if records.len() != bodies.len() {
        return Err("duplicate native water-body id".into());
    }
    let mut surface = Vec::with_capacity(ids.len());
    let mut kinds = Vec::with_capacity(ids.len());
    let mut quantization = 0.0f64;
    let mut over = 0.0f64;
    for (c, (&z, &id)) in elevation.iter().zip(ids).enumerate() {
        if !z.is_finite() || radius + z as f64 <= 0.0 {
            return Err(format!("invalid radial solid height at cell {c}").into());
        }
        if id == u32::MAX {
            surface.push(z);
            kinds.push(0);
            continue;
        }
        let b = records
            .get(&id)
            .ok_or_else(|| format!("unknown water body at cell {c}"))?;
        let stage = finite(b.water_surface_m, "water stage")?;
        let code = match b.kind {
            WaterBodyKindV0::Marine => 1,
            WaterBodyKindV0::Lake => 2,
        };
        let tolerance = 2.0 * f32::EPSILON as f64 * (z as f64).abs().max(stage.abs()).max(1.0);
        if z as f64 - stage > tolerance {
            return Err(format!("native water stage below its solid bed at cell {c}").into());
        }
        if code == 1 && (stage - sea).abs() > tolerance {
            return Err("marine stage differs from declared sea level".into());
        }
        let value = stage as f32;
        if !value.is_finite() || radius + stage <= 0.0 {
            return Err("water stage cannot be represented as a finite radial f32 height".into());
        }
        quantization = quantization.max((value as f64 - stage).abs());
        over = over.max(z as f64 - stage);
        surface.push(value);
        kinds.push(code);
    }
    Ok((surface, kinds, quantization, over))
}
fn top_lithology(s: &StratigraphyDomainV0, n: usize) -> Result<Vec<u8>> {
    s.validate()?;
    if s.nlat as usize * s.nlon as usize != n {
        return Err("stratigraphy shape mismatch".into());
    }
    Ok((0..n)
        .map(|c| {
            (s.layer_offset[c] as usize..s.layer_offset[c + 1] as usize)
                .find(|&k| s.thickness_m[k] > 0.0)
                .map(|k| s.lith[k])
                .unwrap_or(255)
        })
        .collect())
}
fn add_field(
    fields: &mut BTreeMap<String, Value>,
    payloads: &mut Vec<(String, Vec<u8>)>,
    key: &str,
    file: &str,
    ty: &str,
    unit: &str,
    bytes: Vec<u8>,
    source: Value,
) {
    fields.insert(key.into(),json!({"file":file,"type":ty,"unit":unit,"sha256":sha(&bytes),"bytes":bytes.len(),"source":source}));
    payloads.push((file.into(), bytes));
}
fn export(input: &Path, output: &Path) -> Result<Value> {
    if output.exists() {
        return Err(
            "output already exists; choose a new directory to preserve prior results".into(),
        );
    }
    let input = input.canonicalize()?;
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .canonicalize()?;
    let output = parent.join(
        output
            .file_name()
            .ok_or("output must name a new directory")?,
    );
    if output.starts_with(&input) {
        return Err("output must be outside the read-only source world".into());
    }
    let direct = input.join("terrain/artifact.json").is_file();
    let wrapped = input.join("run/terrain/artifact.json").is_file();
    let root = match (direct, wrapped) {
        (true, false) => input.clone(),
        (false, true) => input.join("run"),
        _ => {
            return Err(
                "expected exactly one published world at WORLD_ROOT or WORLD_ROOT/run".into(),
            )
        }
    };
    let terrain_dir = root.join("terrain");
    let hydro_dir = root.join("hydrology");
    let terrain = checked_artifact(&terrain_dir)?;
    let hydro = checked_artifact(&hydro_dir)?;
    if terrain.header.schema != "pf.terrain/v0" || hydro.header.schema != "pf.hydrology/v0" {
        return Err("unsupported native artifact schemas".into());
    }
    require_inputs(&hydro, &[&terrain.digest])?;
    if terrain.header.semantic.seed_path != hydro.header.semantic.seed_path {
        return Err("terrain/hydrology seed mismatch".into());
    }
    let grid = field_grid(&terrain, "elevation_m", Dtype::F32Le, "m")?;
    let (rows, cols, radius) = grid;
    let n = rows * cols;
    // Water body identity is a native per-cell unsigned code; do not infer water from elevation alone.
    let wf = hydro
        .header
        .semantic
        .fields
        .iter()
        .find(|f| f.name == "water_body_id")
        .ok_or("native water-body field missing")?;
    check_grid(
        field_grid(&hydro, "water_body_id", Dtype::U32Le, &wf.units)?,
        grid,
    )?;
    let epoch = number(&terrain.header.domain, "t_my")?;
    let sea = number(&terrain.header.domain, "sea_level_m")?;
    if epoch < 0.0 || number(&hydro.header.domain, "t_my")? != epoch {
        return Err("terrain/hydrology epoch mismatch".into());
    }
    let domain: HydrologyDomainV0 = serde_json::from_value(hydro.header.domain.clone())?;
    domain.validate()?;
    if domain
        .water_bodies
        .iter()
        .any(|b| b.seed_element.0 >= n as u64)
    {
        return Err("water-body seed outside source grid".into());
    }
    let elevation = checked_floats(&terrain, &terrain_dir, "elevation_m", n)?;
    let water_bytes = native_bytes(&hydro, &hydro_dir, "water_body_id", Dtype::U32Le)?;
    let water_ids: Vec<u32> = water_bytes
        .chunks_exact(4)
        .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    let (surface, kinds, quantization, over) =
        water_surface(&elevation, &water_ids, &domain.water_bodies, sea, radius)?;
    let mut lithology = vec![255; n];
    let mut labels = BTreeMap::from([("255".to_string(), "unknown".to_string())]);
    let mut strata_source = Value::Null;
    if let Some(sidecar) = terrain
        .header
        .semantic
        .sidecars
        .iter()
        .find(|s| s.logical_name == "stratigraphy")
    {
        let bytes = std::fs::read(terrain_dir.join(&sidecar.path))?;
        if bytes.len() as u64 != sidecar.byte_length
            || pf_contracts::sim_module::digest_bytes(&bytes) != sidecar.digest
        {
            return Err("stratigraphy sidecar changed after source verification".into());
        }
        let s: StratigraphyDomainV0 = serde_json::from_slice(&bytes)?;
        if s.t_my != epoch || s.nlat as usize != rows || s.nlon as usize != cols {
            return Err("stratigraphy epoch/grid mismatch".into());
        }
        lithology = top_lithology(&s, n)?;
        for (code, kind) in &s.codes.lithology {
            let value = serde_json::to_value(kind)?;
            labels.insert(
                code.to_string(),
                value
                    .as_str()
                    .ok_or("non-string lithology code label")?
                    .into(),
            );
        }
        strata_source = json!({"artifact":"terrain","sidecar":"stratigraphy","sha256":sha(&bytes),"nativeDigest":sidecar.digest});
    }
    let mut fields = BTreeMap::new();
    let mut payloads = Vec::new();
    add_field(
        &mut fields,
        &mut payloads,
        "elevation",
        "elevation.f32",
        "f32le",
        "m",
        f32_bytes(&elevation),
        json!({"artifact":"terrain","field":"elevation_m"}),
    );
    add_field(
        &mut fields,
        &mut payloads,
        "surfaceHeight",
        "surface-height.f32",
        "f32le",
        "m",
        f32_bytes(&surface),
        json!({"artifact":"hydrology","rule":"native body stage on wet cells, original elevation on dry cells"}),
    );
    add_field(
        &mut fields,
        &mut payloads,
        "waterKind",
        "water-kind.u8",
        "u8",
        "code",
        kinds.clone(),
        json!({"artifact":"hydrology","field":"water_body_id","rule":"native body kind; dry sentinel remains dry"}),
    );
    add_field(
        &mut fields,
        &mut payloads,
        "lithology",
        "lithology.u8",
        "u8",
        "code",
        lithology.clone(),
        strata_source.clone(),
    );
    let mut artifacts = BTreeMap::from([
        ("terrain", terrain.digest.clone()),
        ("hydrology", hydro.digest.clone()),
    ]);
    let mut climate_dirs = vec![root.join("climate")];
    if root != input {
        climate_dirs.push(input.join("climate"));
    }
    let present: Vec<_> = climate_dirs.into_iter().filter(|p| p.exists()).collect();
    if present.len() > 1 {
        return Err("ambiguous climate source; retain one coherent climate artifact".into());
    }
    let climate_available = if let Some(dir) = present.first() {
        let climate = checked_artifact(dir)?;
        require_inputs(&climate, &[&terrain.digest, &hydro.digest])?;
        if climate.header.semantic.seed_path != terrain.header.semantic.seed_path {
            return Err("climate seed mismatch".into());
        }
        for (key, name, unit) in [
            ("temperature", "temperature_mean_c", "degC"),
            ("precipitation", "precipitation_mm_per_yr", "mm/yr"),
        ] {
            check_grid(field_grid(&climate, name, Dtype::F32Le, unit)?, grid)?;
            let values = checked_floats(&climate, dir, name, n)?;
            if key == "precipitation" && values.iter().any(|v| *v < 0.0) {
                return Err("negative native precipitation".into());
            }
            add_field(
                &mut fields,
                &mut payloads,
                key,
                &format!("{key}.f32"),
                "f32le",
                unit,
                f32_bytes(&values),
                json!({"artifact":"climate","field":name}),
            );
        }
        artifacts.insert("climate", climate.digest.clone());
        true
    } else {
        false
    };
    let mut limitations=vec!["Precomputed generated data exported from verified native artifacts; exact reproduction of the historical run by the current source is not claimed.".to_string(),
        "Rendering subdivisions do not add geological resolution. The export is not native fine catchment terrain.".into(),
        "sourceResolutionM is north-south spacing; east-west spacing decreases towards the poles.".into(),
        "Lithology identifies the top positive-thickness published coarse layer, not a fine exposed-material or finite-stock qualification.".into(),
        "Native coarse water-body occupancy and stages are retained; subcell shorelines are unresolved.".into(),
        "Continuous scalar fields use cell-centred bilinear interpolation with pole-row means; categorical codes must not be averaged as numeric quantities.".into()];
    limitations.push(if climate_available {"Native climate is available and pins the same terrain and hydrology artifacts."}else{"Climate is unavailable in this source; temperature and precipitation fields are omitted, not synthesized."}.into());
    if strata_source.is_null() {
        limitations.push(
            "No native stratigraphy sidecar is available; lithology is explicitly unknown (255)."
                .into(),
        );
    }
    let mut manifest = json!({"schema":SCHEMA,"title":"Native generated reference world","seed":terrain.header.semantic.seed_path.seed,
        "radiusM":radius,"epochMy":epoch,"seaLevelM":sea,
        "grid":{"rows":rows,"cols":cols,"order":"south-first","sampling":"cell-centres","layout":"row-major","longitude":"east-positive, periodic","latitude":"north-positive"},
        "sourceResolutionM":radius*std::f64::consts::PI/rows as f64,"interpolation":"cell-centred-bilinear-pole-mean",
        "fields":fields,"sourceArtifacts":artifacts,"lithologyCodes":labels,"waterKindCodes":{"0":"dry","1":"marine","2":"lake"},
        "capabilities":{"climate":climate_available,"coarseLithology":!strata_source.is_null(),"nativeFineTerrain":false},
        "summary":{"heightMinM":elevation.iter().copied().fold(f32::INFINITY,f32::min),"heightMaxM":elevation.iter().copied().fold(f32::NEG_INFINITY,f32::max),
            "dryCells":kinds.iter().filter(|&&v|v==0).count(),"marineCells":kinds.iter().filter(|&&v|v==1).count(),"lakeCells":kinds.iter().filter(|&&v|v==2).count(),
            "unknownLithologyCells":lithology.iter().filter(|&&v|v==255).count(),"maxWaterStageF32ErrorM":quantization,"maxPublishedBedAboveStageM":over},
        "limitations":limitations});
    manifest["id"] = format!("sha256:{}", sha(&serde_json::to_vec(&manifest)?)).into();
    std::fs::create_dir(&output)?;
    for (name, bytes) in payloads {
        std::fs::write(output.join(name), bytes)?;
    }
    // Manifest last: an incomplete write is never a complete viewer data source.
    std::fs::write(
        output.join("manifest.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    Ok(manifest)
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 3 {
        return Err(
            "usage: maiforge_viewer_export WORLD_ROOT OUTPUT_DIR (output parent must exist)".into(),
        );
    }
    let m = export(&PathBuf::from(&args[1]), &PathBuf::from(&args[2]))?;
    println!(
        "{}",
        serde_json::to_string_pretty(
            &json!({"id":m["id"],"grid":m["grid"],"fields":m["fields"].as_object().unwrap().keys().collect::<Vec<_>>(),"climate":m["capabilities"]["climate"],"summary":m["summary"]})
        )?
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use pf_contracts::spatial::ElementId;
    fn body(kind: WaterBodyKindV0, stage: f64) -> WaterBodyV0 {
        WaterBodyV0 {
            body_id: 7,
            kind,
            seed_element: ElementId(0),
            water_surface_m: stage,
            area_m2: 1.0,
            outlet_edges: vec![],
        }
    }
    #[test]
    fn water_uses_native_stages_and_preserves_dry_land() {
        let (h, k, q, _) = water_surface(
            &[12.0, 10.0],
            &[u32::MAX, 7],
            &[body(WaterBodyKindV0::Lake, 13.25)],
            0.0,
            6371000.0,
        )
        .unwrap();
        assert_eq!(h, vec![12.0, 13.25]);
        assert_eq!(k, vec![0, 2]);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn invalid_water_and_nonfinite_fields_are_rejected() {
        assert!(water_surface(
            &[1.0],
            &[7],
            &[body(WaterBodyKindV0::Lake, 0.0)],
            0.0,
            6371000.0
        )
        .is_err());
        assert!(water_surface(
            &[0.0],
            &[8],
            &[body(WaterBodyKindV0::Lake, 0.0)],
            0.0,
            6371000.0
        )
        .is_err());
        assert!(water_surface(
            &[0.0],
            &[7],
            &[body(WaterBodyKindV0::Marine, 1.0)],
            0.0,
            6371000.0
        )
        .is_err());
        assert!(water_surface(&[f32::NAN], &[u32::MAX], &[], 0.0, 6371000.0).is_err());
    }
    #[test]
    fn source_references_cannot_escape_an_artifact() {
        for p in ["/absolute", "../outside", "fields/../../outside", ""] {
            assert!(safe_relative(p).is_err(), "{p}");
        }
        assert!(safe_relative("fields/elevation.f32").is_ok());
        assert!(check_grid((2, 4, 10.0), (2, 4, 11.0)).is_err());
    }
}
