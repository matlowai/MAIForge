//! Past-surface climate observations for source origin stories (D-0154).
//! Reuses the same climate, ocean connectivity and routing kernels as the world.
//! Does not replay the world or feed a new forcing into the frozen C7 history.
use crate::{
    climate::{derive_climate, ClimateInputs},
    routing, water,
};
use pf_contracts::{laws::ClimateRecipeV1, terrain::bed_history::ClimateFrame};
use pf_grid::GridSpec;

pub fn reconstruct_frame(
    grid: &GridSpec,
    time_my: f64,
    elevation_m: &[f64],
    sea_level_m: f64,
    recipe: &ClimateRecipeV1,
    elevation_source_digest: String,
) -> Result<ClimateFrame, String> {
    if !time_my.is_finite() || elevation_source_digest.is_empty() {
        return Err("missing historical time/source identity".into());
    }
    let ocean = water::largest_ocean(grid, elevation_m, sea_level_m).map_err(|e| e.to_string())?;
    // Only routing/fill levels are used; no rainfall proxy is reported as climate.
    let route = routing::route(
        grid,
        elevation_m,
        &vec![0.; grid.cells()],
        &ocean,
        sea_level_m,
    )
    .map_err(|e| e.to_string())?;
    let output = derive_climate(
        &ClimateInputs {
            grid,
            elevation_m,
            fill_level_m: &route.fill_level_m,
            sea_level_m,
            ocean: &ocean,
        },
        recipe,
    )
    .map_err(|e| e.to_string())?;
    let f = output.fields;
    Ok(ClimateFrame {
        time_my,
        elevation_source_digest,
        temperature_c: f.temperature_mean_c.iter().map(|&x| x as f32).collect(),
        aridity_index: f.aridity_index.iter().map(|&x| x as f32).collect(),
        elevation_m: elevation_m.iter().map(|&x| x as f32).collect(),
        water_depth_m: route
            .fill_level_m
            .iter()
            .zip(elevation_m)
            .map(|(&w, &z)| (w - z).max(0.) as f32)
            .collect(),
        ocean_distance_m: f
            .upwind_coast_distance_m
            .iter()
            .map(|&x| x as f32)
            .collect(),
        marine: ocean,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn past_surface_changes_coastal_climate_without_changing_recipe() {
        let g = GridSpec::new(18, 36, 6_371_000.).unwrap();
        let mut z = vec![-1000.; g.cells()];
        for row in 0..18 {
            for col in 10..25 {
                z[row * 36 + col] = 100.;
            }
        }
        let recipe = ClimateRecipeV1::coarse_v0_earth();
        let a = reconstruct_frame(&g, 40., &z, 0., &recipe, "source-a".into()).unwrap();
        for row in 0..18 {
            for col in 15..25 {
                z[row * 36 + col] = 2000.;
            }
        }
        let b = reconstruct_frame(&g, 80., &z, 0., &recipe, "source-b".into()).unwrap();
        let i = 9 * 36 + 20;
        assert!(b.temperature_c[i] < a.temperature_c[i] - 5.);
        assert!(a.ocean_distance_m[i] > 0.);
        assert_eq!(a.water_depth_m[9 * 36], 1000.);
        assert_eq!(a.water_depth_m[i], 0.);
    }
}
