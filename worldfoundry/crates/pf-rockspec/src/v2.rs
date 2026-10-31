//! Version 2 reduced transport and material-frame structure. See MODEL_V2.md.
use super::*;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum StructureKind {
    Crystalline,
    Porphyritic,
    Banded,
    Laminated,
    Clastic,
    Vesicular,
    Glass,
    Micaceous,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Structure {
    pub kind: StructureKind,
    pub optical_boundary_efficiency: Estimate,
    pub wavelength_m: Estimate,
    pub contrast: Estimate,
}
impl Structure {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        require(
            positive(self.optical_boundary_efficiency.value)
                && unit(self.optical_boundary_efficiency.value)
                && positive(self.wavelength_m.value)
                && unit(self.contrast.value),
            "invalid v2 structure",
        )
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct StructureSample {
    pub kind: StructureKind,
    pub wavelength_m: f64,
    pub transport_diameter_m: f64,
    /// Homogenized response at this footprint, before coating and wetness.
    pub mean_reflectance_linear: [f64; 3],
    pub escape_probability: f64,
    pub unresolved_slope_variance: f64,
    pub contrast_filter: f64,
    /// Bounded, resolved optical normal. Not collision/displacement geometry.
    pub normal_material: [f64; 3],
    pub flake_normal_material: [f64; 3],
    pub flake_specular_fraction: f64,
    pub qualification: String,
}
fn hash(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
fn uniform(x: u64) -> f64 {
    (hash(x) >> 11) as f64 / ((1_u64 << 53) as f64)
}
fn seed(m: &Material<'_>) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    // Length-delimited stable material address; NEVER tile, surface, chart or camera.
    for s in [&m.reference.body_id, &m.reference.realization_id] {
        for b in (s.len() as u64).to_le_bytes().into_iter().chain(s.bytes()) {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
    }
    h
}
fn cell_seed(base: u64, cell: [i64; 3]) -> u64 {
    hash(
        base ^ hash(cell[0] as u64)
            ^ hash((cell[1] as u64).wrapping_add(0x12345678))
            ^ hash((cell[2] as u64).wrapping_add(0x76543210)),
    )
}
fn rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [w, x, y, z] = q;
    [
        (1. - 2. * (y * y + z * z)) * v[0]
            + 2. * (x * y - w * z) * v[1]
            + 2. * (x * z + w * y) * v[2],
        2. * (x * y + w * z) * v[0]
            + (1. - 2. * (x * x + z * z)) * v[1]
            + 2. * (y * z - w * x) * v[2],
        2. * (x * z - w * y) * v[0]
            + 2. * (y * z + w * x) * v[1]
            + (1. - 2. * (x * x + y * y)) * v[2],
    ]
}
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt();
    v.map(|x| x / n)
}
fn add_normal(n: [f64; 3], perturb: [f64; 3], strength: f64) -> [f64; 3] {
    let along = dot(n, perturb);
    normalize(std::array::from_fn(|i| {
        n[i] + strength * (perturb[i] - along * n[i])
    }))
}
/// Cavity multiple-return closure. g is escape probability per diffuse return.
/// Includes wall interreflection. g=1 is flat and R=0/1 are exact limits.
pub fn cavity_reflectance(r: f64, g: f64) -> f64 {
    r * g / (1. - (1. - g) * r)
}

pub(crate) fn prepare(
    m: &Material<'_>,
    s: &SurfaceState,
    f: f64,
    optical_grain: f64,
    mix: &mut Mixture,
) -> Result<StructureSample, Error> {
    let desc = m
        .spec
        .structure_v2
        .as_ref()
        .ok_or_else(|| Error::Invalid("missing v2 structure".into()))?;
    desc.validate()?;
    let kind = desc.kind;
    let q = m.bulk.fabric_to_material;
    let p = rotate([q[0], -q[1], -q[2], -q[3]], s.material_coordinate_m);
    let n = rotate([q[0], -q[1], -q[2], -q[3]], s.normal_material);
    let d = desc.wavelength_m.value;
    // Joint relief uses the same fracture slope, spacing, persistence and plane normals
    // as the material descriptor. No per-rock scale/albedo coefficients.
    let mut variance = 0.;
    for j in &m.spec.joint_sets {
        let polar = j.normal_polar_deg.value.to_radians();
        let az = j.normal_azimuth_deg.value.to_radians();
        let jn = [polar.sin() * az.cos(), polar.sin() * az.sin(), polar.cos()];
        let crossing = 1. - dot(n, jn).powi(2);
        let spacing = j.spacing_m.value;
        let support = j.persistence_m.value / (j.persistence_m.value + spacing);
        let relief = 0.5
            * (m.catalogue.model.fracture_slope.value * (1. - mix.glass)
                + m.catalogue.model.glass_slope.value * mix.glass);
        variance += crossing * support * relief.powi(2) * f * f / (f * f + spacing * spacing);
        variance +=
            crossing * (j.aperture_m.value / spacing).powi(2) * f * f / (f * f + spacing * spacing);
    }
    let aspect = m.spec.grain.aspect_ratio.value;
    let grain_slope = 0.18 * (aspect - 1.) / aspect + 0.4 * m.bulk.interparticle_porosity;
    let relief_scale = m.spec.grain.d90_m.value;
    variance += grain_slope.powi(2) * f * f / (f * f + relief_scale * relief_scale);
    let ves = m.spec.vesicularity.value.min(m.bulk.intrinsic_porosity);
    // Deep vesicle escape: quarter opening solid angle (declared cavity approximation).
    let cavity_g = 1. - 0.75 * ves;
    let g = cavity_g / (1. + variance);
    mix.reflectance = mix.reflectance.map(|r| cavity_reflectance(r, g));
    // Wall normals in vesicular material are broad even when glass itself is smooth.
    variance += ves * 0.65 + 0.35 * (m.bulk.intrinsic_porosity - ves);
    let mean = mix.reflectance;
    // Compact filter: below-grid structure is exactly homogenized after 2 wavelengths.
    let filter = (1. - (f / (2. * d)).powi(2)).max(0.).powi(2);
    let mut normal = s.normal_material;
    let mut flake_normal = s.normal_material;
    let mut flake_fraction = 0.;
    if filter > 0. && kind != StructureKind::Glass {
        require(
            p.iter().all(|x| (x / d).abs() < 9e14),
            "material coordinate exceeds structure precision domain",
        )?;
        let base = seed(m);
        let scaled = p.map(|x| x / d);
        let cell = scaled.map(|x| x.floor() as i64);
        let mut best = f64::INFINITY;
        let mut chosen = 0_u64;
        let mut offset = [0.; 3];
        // Jittered material-space Voronoi crystals/clasts. Broad clast population
        // is represented by descriptor d90 support; NOT generated collectible grains.
        for z in -1..=1 {
            for y in -1..=1 {
                for x in -1..=1 {
                    let key = cell_seed(base, [cell[0] + x, cell[1] + y, cell[2] + z]);
                    let jitter = [
                        uniform(key),
                        uniform(key.wrapping_add(1)),
                        uniform(key.wrapping_add(2)),
                    ];
                    let v = [
                        scaled[0] - (cell[0] + x) as f64 - jitter[0],
                        scaled[1] - (cell[1] + y) as f64 - jitter[1],
                        scaled[2] - (cell[2] + z) as f64 - jitter[2],
                    ];
                    let spread = (m.spec.grain.d90_m.value / m.spec.grain.d10_m.value).ln() / 2.563;
                    let size = (spread * (uniform(key.wrapping_add(9)) - 0.5))
                        .exp()
                        .clamp(0.5, 1.5);
                    let dist = dot(v, v) / (size * size);
                    if dist < best {
                        best = dist;
                        chosen = key;
                        offset = v;
                    }
                }
            }
        }
        let u = if matches!(kind, StructureKind::Banded | StructureKind::Laminated) {
            (p[2] / d + uniform(base) + 0.055 * (p[0] / (d * 3.)).sin() + 0.025 * (p[1] / d).sin())
                .rem_euclid(1.)
        } else {
            uniform(chosen.wrapping_add(3))
        };
        let total: f64 = m
            .bulk
            .mineral_mode
            .iter()
            .map(|v| {
                v.mass_fraction
                    / m.catalogue
                        .mineral(&v.mineral)
                        .expect("validated modes")
                        .density_kg_m3
                        .value
            })
            .sum();
        let mut phases = Vec::new();
        let mut areal = [0.; 3];
        for part in &m.bulk.mineral_mode {
            let mineral = m.catalogue.mineral(&part.mineral)?;
            let weight = part.mass_fraction / mineral.density_kg_m3.value / total;
            let single = mixture(
                m.catalogue,
                &[ConstituentFraction {
                    mineral: part.mineral.clone(),
                    mass_fraction: 1.,
                }],
                optical_grain,
            )?;
            let colour = single.reflectance.map(|r| cavity_reflectance(r, g));
            for j in 0..3 {
                areal[j] += weight * colour[j];
            }
            phases.push((
                weight,
                colour,
                mineral.lustre == "pearly",
                part.mineral.as_str(),
            ));
        }
        // Sort mineral response for coherent dark/light segregation in laminae.
        if matches!(kind, StructureKind::Banded | StructureKind::Laminated) {
            phases.sort_by(|a, b| a.1[1].total_cmp(&b.1[1]));
        }
        let mut cumulative = 0.;
        let mut picked = phases
            .last()
            .ok_or_else(|| Error::Invalid("empty phases".into()))?;
        for phase in &phases {
            cumulative += phase.0;
            if u < cumulative {
                picked = phase;
                break;
            }
        }
        let mut amplitude = desc.contrast.value;
        // Upper-decile coarse population for porphyritic specimens; remaining
        // groundmass is homogenized at this population scale.
        let coarse = kind != StructureKind::Porphyritic || uniform(chosen.wrapping_add(10)) < 0.1;
        if !coarse {
            amplitude = 0.;
        }
        for phase in &phases {
            for j in 0..3 {
                let delta = phase.1[j] - areal[j];
                if delta < 0. {
                    amplitude = amplitude.min(mean[j] / (-delta));
                }
                if delta > 0. {
                    amplitude = amplitude.min((1. - mean[j]) / delta);
                }
            }
        }
        // Zero-mean optical contrast around the KM bulk response. No mineral or ore
        // is admitted by this diagnostic projection. Exact far-field value is mean.
        for j in 0..3 {
            mix.reflectance[j] =
                (mean[j] + amplitude * filter * (picked.1[j] - areal[j])).clamp(0., 1.);
        }
        let slope = match kind {
            StructureKind::Clastic => 0.5,
            StructureKind::Vesicular => 0.8,
            StructureKind::Micaceous => 0.08,
            _ => 0.16,
        };
        let mut delta = offset;
        if matches!(kind, StructureKind::Banded | StructureKind::Laminated) {
            delta = [0., 0., (u * std::f64::consts::TAU).cos()];
        }
        if kind == StructureKind::Vesicular {
            // Random compact spherical pores wholly contained in material cells.
            // Rotate the isotropic pore lattice to avoid aligned planar sections.
            let pq = [0.7, 0.3, 0.5, 0.2].map(|x| x / 0.87_f64.sqrt());
            let pp = rotate(pq, scaled);
            let pc = pp.map(|x| x.floor() as i64);
            let pk = cell_seed(base.wrapping_add(71), pc);
            let radius = 0.2 + 0.2 * uniform(pk.wrapping_add(3));
            let pv: [f64; 3] = std::array::from_fn(|i| {
                pp[i] - pc[i] as f64 - (0.4 + 0.2 * uniform(pk.wrapping_add(i as u64)))
            });
            let t = (1. - dot(pv, pv) / (radius * radius)).max(0.);
            let pore = t * t;
            // Exact volume integral of (1-r²/a²)² over a sphere, averaged over
            // uniform radius [0.2,0.4]. A planar finite patch has sampling error.
            let pore_mean = 32. * std::f64::consts::PI / 105. * (0.4_f64.powi(4) - 0.2_f64.powi(4))
                / (4. * 0.2);
            for (j, mean_channel) in mean.iter().enumerate() {
                let amplitude = (mean_channel / (1. - pore_mean))
                    .min((1. - mean_channel) / pore_mean)
                    * 0.9
                    * ves
                    * filter;
                mix.reflectance[j] = mean_channel + amplitude * (pore_mean - pore);
            }
            let gradient = pv.map(|x| -4. * t * x / (radius * radius));
            delta = rotate([pq[0], -pq[1], -pq[2], -pq[3]], gradient).map(|x| x * 0.18);
        }
        normal = add_normal(
            s.normal_material,
            rotate(q, delta),
            slope * filter * if coarse { 1. } else { 0.1 },
        );
        if picked.2
            && matches!(
                kind,
                StructureKind::Micaceous | StructureKind::Laminated | StructureKind::Crystalline
            )
        {
            let flake = normalize([
                0.9 * (uniform(chosen.wrapping_add(4)) - 0.5),
                0.9 * (uniform(chosen.wrapping_add(5)) - 0.5),
                1.,
            ]);
            flake_normal = rotate(q, flake);
            flake_fraction = m.catalogue.model.sparkle_fraction.value * filter;
        }
    }
    Ok(StructureSample {kind,wavelength_m:d,transport_diameter_m:optical_grain,mean_reflectance_linear:mean,escape_probability:g,unresolved_slope_variance:variance,contrast_filter:filter,normal_material:normal,flake_normal_material:flake_normal,flake_specular_fraction:flake_fraction,qualification:"conditional optical structure in rest-material frame; not stock/geometry; approximate compact filter; joint closure is not a full Hapke BRDF".into()})
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ExposureScenario {
    WetTemperate,
    AridVarnish,
}
/// Explicit manufactured exposure state for diagnostics, NOT a weathering rate law.
/// Arid varnish includes externally supplied Mn/Fe dust. Callers must admit that
/// forcing/deposition separately in production. Never called implicitly by sampling.
pub fn diagnostic_weathering(
    m: &Material<'_>,
    mut s: SurfaceState,
    scenario: ExposureScenario,
) -> Result<SurfaceState, Error> {
    require(
        m.catalogue.schema_version == "pf-rockspec/2",
        "v2 weathering requires v2 catalogue",
    )?;
    let family = m.spec.weathering_family.as_str();
    s.mechanism = FaceMechanism::Inherited;
    s.exposure_elapsed_s = None;
    s.forcing_record=format!("manufactured {scenario:?} exposure for {family}; supplied residue/dust; no age/rate or stock transfer claimed");
    s.rms_slope_history_rad = 0.08;
    s.abrasion = 0.05;
    let (mode, coverage, thickness): (Vec<(&str, f64)>, f64, f64) = if family == "organic_oxidation"
    {
        (
            vec![("organic_matter", 0.97), ("goethite", 0.03)],
            0.25,
            0.00004,
        )
    } else if scenario == ExposureScenario::AridVarnish {
        (
            vec![("pyrolusite", 0.45), ("goethite", 0.30), ("illite", 0.25)],
            0.55,
            0.00006,
        )
    } else {
        match family {
            "mafic_oxidation" | "serpentinization" => (
                vec![("goethite", 0.6), ("hematite", 0.15), ("smectite", 0.25)],
                0.4,
                0.00005,
            ),
            "carbonate_solution" => (
                vec![("calcite", 0.55), ("illite", 0.4), ("organic_matter", 0.05)],
                0.25,
                0.00006,
            ),
            "organic_oxidation" => unreachable!(),
            "evaporite_solution" => (vec![("illite", 1.)], 0., 0.),
            "slaking" => (
                vec![
                    ("illite", 0.65),
                    ("goethite", 0.15),
                    ("organic_matter", 0.2),
                ],
                0.25,
                0.00005,
            ),
            "silica_residual" => (
                vec![
                    ("quartz", 0.65),
                    ("illite", 0.22),
                    ("goethite", 0.1),
                    ("organic_matter", 0.03),
                ],
                0.2,
                0.00005,
            ),
            "silicate_hydrolysis"
            | "grus_disaggregation"
            | "residual_weathering"
            | "granular_disaggregation"
            | "glass_hydration" => (
                vec![("kaolinite", 0.55), ("illite", 0.3), ("goethite", 0.15)],
                0.18,
                0.00005,
            ),
            _ => {
                return Err(Error::Invalid(format!(
                    "unsupported weathering family {family}"
                )))
            }
        }
    };
    s.coating = if coverage == 0. {
        None
    } else {
        Some(CoatingLayer {
            source_event: format!("diagnostic-v2-{scenario:?}-{family}"),
            mineral_mode: mode
                .into_iter()
                .map(|(k, v)| ConstituentFraction {
                    mineral: k.into(),
                    mass_fraction: v,
                })
                .collect(),
            grain_diameter_m: 0.0002,
            thickness_m: thickness,
            coverage,
            extinction_m_inv: m.catalogue.model.rind_extinction_m_inv.value,
        })
    };
    Ok(s)
}
