//! D-0139 / D-0154. Cited rock-family relationships with deliberately coarse
//! global P-T boundaries. An observation over retained C7 protolith/burial;
//! NEVER mutates native stacks, stock quantities, mineral modes or lineage.
use crate::stratigraphy::*;
use pf_contracts::terrain::bed_history::MetamorphicFacts;

/// ColdSubduction requires evidence of this material descending in a slab.
/// An arc source, nearby trench or surface boundary flag is NOT that evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalRegime {
    Regional,
    ColdSubduction,
}

pub const LAW: &str = "protolith-burial-pt-reconstruction-v1";
/// Pressure is lithostatic rho*g*z. Gradients and boundaries are declared
/// world-wide approximations, not thermobarometry. No guessed pressure from grade.
pub fn reconstruct(
    protolith: u8,
    burial_m: Option<f64>,
    gravity: f64,
    crust_density: f64,
    regime: ThermalRegime,
) -> Result<MetamorphicFacts, String> {
    if protolith as usize >= NLITH
        || !gravity.is_finite()
        || gravity <= 0.
        || !crust_density.is_finite()
        || crust_density <= 0.
        || burial_m.is_some_and(|v| !v.is_finite() || v < 0.)
    {
        return Err("invalid metamorphic protolith/depth/planet constants".into());
    }
    let cold = regime == ThermalRegime::ColdSubduction;
    let mut r = MetamorphicFacts {
        law: LAW.into(), protolith, peak_burial_m: burial_m,
        pressure_gpa: None, temperature_c: None, temperature_range_c: None,
        regime: if cold {"cold_subduction"} else {"regional_burial"}.into(),
        grade: "unknown".into(), rock_name: None, assemblage: vec![],
        status: "unknown: no retained contribution burial history".into(),
        assumptions: vec![
            "C7-derived reconstruction; no native P-T-t or fluid history".into(),
            "P=rho*g*z; source maximum burial is layer-bottom depth, an upper-depth proxy for the whole contribution".into(),
            "T=15 C + gradient*depth; regional 25 (20..30) C/km, cold slab 10 (7..15) C/km; global assumptions, not measured confidence intervals".into(),
            "Qualitative assemblage only; no phase proportions, devolatilization, melt extraction or inventory conversion. Existing mass, elements and parent lineage are retained unchanged".into(),
            "Texture and reaction progress depend on strain, time, fluids and chemistry unrecorded by C7; names are conditional, not measured".into(),
        ],
        citations: vec![
            "https://www.maine.gov/DACF/mgs/education/lessons/act25.pdf".into(),
            "https://apps.usgs.gov/thesaurus/thesaurus-full.php?thcode=4".into(),
            "https://doi.org/10.1016/0040-1951(73)90006-1".into(),
            "https://doi.org/10.1046/j.1440-1738.1998.00180.x".into(),
        ],
    };
    let Some(z) = burial_m else {
        return Ok(r);
    };
    let p = crust_density * gravity * z / 1e9;
    let t = 15. + z / 1000. * if cold { 10. } else { 25. };
    r.pressure_gpa = Some(p);
    r.temperature_c = Some(t);
    r.temperature_range_c = Some(if cold {
        [15. + 7. * z / 1000., 15. + 15. * z / 1000.]
    } else {
        [15. + 20. * z / 1000., 15. + 30. * z / 1000.]
    });
    if t > 800. || p > 3. {
        r.status =
            "outside solid-state approximation: melting/ultrahigh-pressure evolution unresolved"
                .into();
        return Ok(r);
    }
    if protolith == LITH_BASEMENT || protolith == LITH_METAMORPHIC {
        r.status =
            "unknown: original protolith unresolved; generic metamorphic code is not a protolith"
                .into();
        return Ok(r);
    }
    let mafic = protolith == LITH_BASALT || protolith == LITH_OCEANIC;
    r.grade = if cold && p >= 1.5 && t >= 450. {
        "high_pressure_eclogite"
    } else if cold && p >= 0.6 && (200. ..=550.).contains(&t) {
        "high_pressure_blueschist"
    } else if t < 200. {
        "diagenetic_or_unmetamorphosed"
    } else if t < 350. {
        "very_low"
    } else if t < 450. {
        "low"
    } else if t < 650. {
        "medium"
    } else {
        "high"
    }
    .into();
    let (rock, phases): (&str, &[&str]) = if t < 200. {
        (LITH_NAMES[protolith as usize], &[])
    } else if mafic && r.grade == "high_pressure_eclogite" {
        ("eclogite", &["garnet", "omphacite"])
    } else if mafic && r.grade == "high_pressure_blueschist" {
        (
            "blueschist",
            &["glaucophane", "lawsonite_or_epidote", "chlorite"],
        )
    } else {
        match protolith {
            LITH_SHALE if t < 350. => ("slate", &["quartz", "white_mica", "chlorite"]),
            LITH_SHALE if t < 450. => {
                ("phyllite", &["quartz", "white_mica", "chlorite_or_biotite"])
            }
            LITH_SHALE if t < 650. => ("schist", &["quartz", "mica", "garnet_possible"]),
            LITH_SHALE => (
                "gneiss",
                &["quartz", "feldspar", "biotite", "aluminosilicate_possible"],
            ),
            LITH_CARBONATE => ("marble", &["calcite_or_dolomite"]),
            LITH_SANDSTONE => ("quartzite", &["quartz"]),
            LITH_BASALT | LITH_OCEANIC if t < 300. => (
                "metabasalt",
                &["chlorite", "albite", "prehnite_or_pumpellyite_possible"],
            ),
            LITH_BASALT | LITH_OCEANIC if t < 500. => (
                "greenschist",
                &["chlorite", "actinolite", "epidote", "albite"],
            ),
            LITH_BASALT | LITH_OCEANIC => ("amphibolite", &["hornblende", "plagioclase"]),
            LITH_GRANITE if t >= 550. => ("gneiss", &["quartz", "feldspar", "mica"]),
            LITH_GRANITE => ("metagranite", &["quartz", "feldspar", "mica"]),
            LITH_CONGLOMERATE => ("metaconglomerate", &[]),
            _ => {
                r.status = "unsupported bulk composition for metamorphic assemblage".into();
                return Ok(r);
            }
        }
    };
    r.rock_name = Some(rock.into());
    r.assemblage = phases.iter().map(|x| (*x).into()).collect();
    r.status = "conditional protolith-aware P-T approximation".into();
    if protolith == LITH_SANDSTONE {
        r.assumptions.push("Quartzite presumes quartz-rich sandstone; feldspathic/muddy compositions need a specific assay".into());
    }
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn protolith_and_slab_evidence_change_the_story() {
        let m = |l, z, r| reconstruct(l, Some(z), 9.81, 2800., r).unwrap();
        let regional = ThermalRegime::Regional;
        assert_eq!(
            m(LITH_SHALE, 12000., regional).rock_name.as_deref(),
            Some("slate")
        );
        assert_eq!(
            m(LITH_SHALE, 16000., regional).rock_name.as_deref(),
            Some("phyllite")
        );
        assert_eq!(
            m(LITH_SHALE, 20000., regional).rock_name.as_deref(),
            Some("schist")
        );
        assert_eq!(
            m(LITH_SHALE, 27000., regional).rock_name.as_deref(),
            Some("gneiss")
        );
        assert_eq!(
            m(LITH_CARBONATE, 20000., regional).rock_name.as_deref(),
            Some("marble")
        );
        assert_eq!(
            m(LITH_SANDSTONE, 20000., regional).rock_name.as_deref(),
            Some("quartzite")
        );
        assert_eq!(
            m(LITH_BASALT, 14000., regional).rock_name.as_deref(),
            Some("greenschist")
        );
        assert_eq!(
            m(LITH_BASALT, 22000., regional).rock_name.as_deref(),
            Some("amphibolite")
        );
        assert_eq!(
            m(LITH_BASALT, 30000., ThermalRegime::ColdSubduction)
                .rock_name
                .as_deref(),
            Some("blueschist")
        );
        assert_eq!(
            m(LITH_BASALT, 60000., ThermalRegime::ColdSubduction)
                .rock_name
                .as_deref(),
            Some("eclogite")
        );
        assert_ne!(
            m(LITH_BASALT, 30000., regional).rock_name.as_deref(),
            Some("blueschist")
        );
        assert!(m(LITH_METAMORPHIC, 20000., regional).rock_name.is_none());
        assert!(reconstruct(LITH_SHALE, None, 9.81, 2800., regional)
            .unwrap()
            .rock_name
            .is_none());
        assert!(m(LITH_SHALE, 60000., regional).rock_name.is_none());
        assert!(reconstruct(LITH_SHALE, Some(-1.), 9.81, 2800., regional).is_err());
    }
}
