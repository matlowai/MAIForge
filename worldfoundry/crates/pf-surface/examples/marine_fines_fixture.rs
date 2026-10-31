//! Q25 controlled native surface-step fixture, not a C7 replay or new world.
//! Usage: marine_fines_fixture C7_ROOT OUTPUT_JSON
use pf_grid::GridSpec;
use pf_surface::erosion::marine_fines::MarineFinesParams;
use pf_surface::erosion::retain::{
    erode_retain_step_with_lakes, erode_retain_step_with_marine_fines, RetentionParams,
};
use pf_surface::erosion::{DiffusionStencil, ErosionParams, ErosionStep};
use pf_surface::routing::route;
use pf_surface::sediment::CompensatedColumn;
use serde_json::{json, Value};
use std::error::Error;
use std::io::Write;
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const SOURCE_ROWS: usize = 288;
const SOURCE_COLS: usize = 576;
const BLOCK: usize = 4;

fn read_f32(path: &Path) -> Result<Vec<f64>> {
    let bytes = std::fs::read(path)?;
    if bytes.len() != SOURCE_ROWS * SOURCE_COLS * 4 {
        return Err(format!("{}: expected 576x288 f32", path.display()).into());
    }
    let values: Vec<f64> = bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64)
        .collect();
    if values.iter().any(|v| !v.is_finite()) {
        return Err(format!("{}: nonfinite field", path.display()).into());
    }
    Ok(values)
}

fn average_blocks(values: &[f64], source: &GridSpec) -> Result<Vec<f64>> {
    let mut result = Vec::with_capacity(SOURCE_ROWS * SOURCE_COLS / (BLOCK * BLOCK));
    for i in 0..SOURCE_ROWS / BLOCK {
        for j in 0..SOURCE_COLS / BLOCK {
            let mut sum = 0.0;
            let mut area = 0.0;
            for di in 0..BLOCK {
                let row = i * BLOCK + di;
                let a = source.cell_area_m2(row)?;
                for dj in 0..BLOCK {
                    sum += a * values[row * SOURCE_COLS + j * BLOCK + dj];
                    area += a;
                }
            }
            result.push(sum / area);
        }
    }
    Ok(result)
}

fn required_number(value: &Value, pointer: &str) -> Result<f64> {
    value
        .pointer(pointer)
        .and_then(Value::as_f64)
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("missing finite recipe field {pointer}").into())
}

fn metrics(
    step: &ErosionStep,
    initial_sediment: &[f64],
    elevation: &[f64],
    areas: &[f64],
) -> Result<Value> {
    let mut deep_area = 0.0;
    let mut bare_area = 0.0;
    let mut volume = 0.0;
    let mut added_volume = 0.0;
    for c in 0..areas.len() {
        let after = initial_sediment[c] + step.sediment_added_m[c] - step.sediment_removed_m[c];
        if !after.is_finite() || after < -1e-7 {
            return Err(format!("invalid final sediment at {c}: {after}").into());
        }
        if elevation[c] < -2500.0 {
            deep_area += areas[c];
            if after == 0.0 {
                bare_area += areas[c];
            }
            volume += after * areas[c];
            added_volume += step.sediment_added_m[c] * areas[c];
        }
    }
    if deep_area <= 0.0 {
        return Err("fixture contains no depth<-2500m receivers".into());
    }
    let residual = step.eroded_m3 - step.deposited_m3 - step.exported_m3;
    let array_added: f64 = step
        .sediment_added_m
        .iter()
        .zip(areas)
        .map(|(h, a)| h * a)
        .sum();
    let array_residual = array_added - step.eroded_m3;
    let tolerance = 1.0 + 1e-8 * step.eroded_m3;
    if residual.abs() > tolerance || array_residual.abs() > tolerance {
        return Err(format!("native mass failure: {residual}, {array_residual}").into());
    }
    Ok(
        json!({"deep_mask":"initial elevation < -2500 m, all crust kinds",
        "deep_area_m2":deep_area,"deep_bare_area_fraction":bare_area/deep_area,
        "deep_mean_sediment_m":volume/deep_area,"deep_mean_added_m":added_volume/deep_area,
        "eroded_m3":step.eroded_m3,"land_lake_deposited_m3":step.deposited_m3,
        "marine_received_m3":step.exported_m3,"closure_residual_m3":residual,
        "closure_relative":residual.abs()/step.eroded_m3.max(1.0),
        "added_array_residual_m3":array_residual}),
    )
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: marine_fines_fixture C7_ROOT OUTPUT_JSON".into());
    }
    let world = Path::new(&args[1]);
    let out = Path::new(&args[2]);
    let recipe: Value = serde_json::from_slice(&std::fs::read(world.join("recipe.json"))?)?;
    let profile = &recipe["inputs"]["profile"];
    let radius = required_number(profile, "/planet/radius_m")?;
    let source = GridSpec::new(SOURCE_ROWS, SOURCE_COLS, radius)?;
    let grid = GridSpec::new(SOURCE_ROWS / BLOCK, SOURCE_COLS / BLOCK, radius)?;
    let fields = world.join("terrain/fields");
    let elevation = average_blocks(&read_f32(&fields.join("elevation_m.f32"))?, &source)?;
    let sediment = average_blocks(
        &read_f32(&fields.join("sediment_thickness_m.f32"))?,
        &source,
    )?;
    let available = average_blocks(&read_f32(&fields.join("crust_thickness_m.f32"))?, &source)?;
    let kind = std::fs::read(fields.join("crust_kind.u8"))?;
    if kind.len() != SOURCE_ROWS * SOURCE_COLS || kind.iter().any(|&v| v > 1) {
        return Err("invalid kind field, expected 576x288 u8 oceanic=0/continental=1".into());
    }
    let continental_fraction =
        average_blocks(&kind.iter().map(|&v| v as f64).collect::<Vec<_>>(), &source)?;
    let continental: Vec<bool> = continental_fraction.iter().map(|&f| f >= 0.5).collect();
    // Published elevations are relative to sea level; resampling changes coasts.
    let ocean: Vec<bool> = elevation.iter().map(|&z| z < 0.0).collect();
    let n = grid.cells();
    let stencil = DiffusionStencil::build(&grid)?;
    let runoff = required_number(profile, "/surface/uniform_runoff_m_per_yr")? * 1e6;
    let routing = route(&grid, &elevation, &vec![runoff; n], &ocean, 0.0)?;
    let column = CompensatedColumn {
        mantle_density_kg_m3: required_number(profile, "/planet/mantle_density_kg_m3")?,
        sediment_density_kg_m3: required_number(profile, "/planet/crust_density_kg_m3")?,
        water_density_kg_m3: required_number(
            profile,
            "/two_column_isostasy/seawater_density_kg_m3",
        )?,
    };
    let retention_recipe = &profile["orogeny"]["sediment_retention"];
    let mut retention = RetentionParams::v1(column);
    retention.base_level_relief_scale_m =
        required_number(retention_recipe, "/base_level_relief_scale_m")?;
    retention.taper_exponent = required_number(retention_recipe, "/base_level_taper_exponent")?;
    retention.alluvial_grade_slope = required_number(retention_recipe, "/alluvial_grade_slope")?;
    retention.shelf_wave_base_depth_m =
        required_number(retention_recipe, "/shelf_wave_base_depth_m")?;
    let erosion = ErosionParams {
        k_per_my: required_number(
            profile,
            "/pinned_universe/constants/stream_power_erodibility_per_yr/value",
        )? * 1e6,
        m: 0.5,
        n: 1.0,
        diffusivity_m2_per_my: required_number(
            profile,
            "/pinned_universe/constants/hillslope_diffusivity_m2_per_yr/value",
        )? * 1e6,
        stability_limit: 0.5,
    };
    let fines = MarineFinesParams {
        fraction: 0.25,
        min_depth_m: 2500.0,
    };
    let dt = 0.1;
    let start = std::time::Instant::now();
    let off = erode_retain_step_with_lakes(
        &grid,
        &stencil,
        &routing,
        &elevation,
        &sediment,
        &ocean,
        0.0,
        dt,
        &erosion,
        Some(&available),
        &continental,
        &retention,
        None,
        false,
        None,
    )?;
    let off_seconds = start.elapsed().as_secs_f64();
    let start = std::time::Instant::now();
    let on = erode_retain_step_with_marine_fines(
        &grid,
        &stencil,
        &routing,
        &elevation,
        &sediment,
        &ocean,
        0.0,
        dt,
        &erosion,
        Some(&available),
        &continental,
        &retention,
        None,
        false,
        None,
        Some(&fines),
    )?;
    let on_seconds = start.elapsed().as_secs_f64();
    for c in 0..n {
        if on.0.bedrock_removed_m[c] > available[c] + 1e-7
            || off.0.bedrock_removed_m[c] > available[c] + 1e-7
        {
            return Err(format!("finite bedrock overspent at {c}").into());
        }
    }
    let off_metrics = metrics(&off.0, &sediment, &elevation, &stencil.areas_m2)?;
    let on_metrics = metrics(&on.0, &sediment, &elevation, &stencil.areas_m2)?;
    let result = json!({
        "schema":"vz.surface-test.marine-fines/v1",
        "kind":"controlled one-step native surface fixture, not C7 replay or world run",
        "world_input":world,"shape":[grid.nlat(),grid.nlon()],"source_shape":[SOURCE_ROWS,SOURCE_COLS],
        "resampling":"4x4 area-weighted means, south-first; continental area fraction >=0.5; ocean from mean elevation<0",
        "limitations":["Endpoint fields, not historical replay","No tectonic advance, isostatic re-solve, climate update, shelf grading or stratigraphy capture","Native base retention step uses n=1; C7's nonlinear/subgrid and incision-threshold preparation are not reproduced","No pelagic mass or full grain classes","No production recipe activation"],
        "parameters":{"dt_my":dt,"k_per_my":erosion.k_per_my,"m":erosion.m,"n":erosion.n,
            "diffusivity_m2_per_my":erosion.diffusivity_m2_per_my,"stability_limit":erosion.stability_limit,
            "fines":fines,"retention_recipe":retention_recipe,"column":{"mantle":column.mantle_density_kg_m3,"sediment":column.sediment_density_kg_m3,"water":column.water_density_kg_m3}},
        "initial_elevation_m":elevation,"initial_sediment_m":sediment,"initial_bedrock_m":available,
        "ocean":ocean,"continental":continental,"cell_area_m2":stencil.areas_m2,
        "off":{"metrics":off_metrics,"seconds":off_seconds,"sediment_added_m":off.0.sediment_added_m,
            "sediment_removed_m":off.0.sediment_removed_m,"bedrock_removed_m":off.0.bedrock_removed_m},
        "on":{"metrics":on_metrics,"seconds":on_seconds,"sediment_added_m":on.0.sediment_added_m,
            "sediment_removed_m":on.0.sediment_removed_m,"bedrock_removed_m":on.0.bedrock_removed_m,"fines_report":on.3}
    });
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let file = std::fs::File::create(out)?;
    let mut writer = std::io::BufWriter::new(file);
    serde_json::to_writer(&mut writer, &result)?;
    writer.flush()?;
    println!(
        "{}",
        serde_json::to_string_pretty(
            &json!({"output":out,"off":result["off"]["metrics"],"on":result["on"]["metrics"],"seconds":{"off":off_seconds,"on":on_seconds}})
        )?
    );
    Ok(())
}
