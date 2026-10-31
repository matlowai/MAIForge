//! Passive mixture transport on the actual stratigraphy mutation graph.
//! Compression merges consecutive formation intervals of the same carrier/process.
//! Erosion and extension retain homogeneous fractions within compacted layers;
//! basement reconciliation is UNKNOWN, never attributed to adjacent magmatism.
use crate::{
    layer_evolution::{LayerChange, BASEMENT},
    stratigraphy::*,
};
pub use pf_contracts::layer_provenance::{Contribution, LayerProvenance, Setting, StackProvenance};

fn single(
    h: f64,
    setting: Setting,
    carrier: u64,
    process: &str,
    lith: u8,
    t0: f64,
    t1: f64,
) -> LayerProvenance {
    LayerProvenance {
        contributions: if h > 0. {
            vec![Contribution {
                thickness_m: h,
                setting,
                episode: format!("carrier:{carrier}:{process}:{:016x}", t0.to_bits()),
                start_my: t0,
                end_my: t1,
                protolith: lith,
                peak_grade: None,
                peak_burial_m: None,
            }]
        } else {
            vec![]
        },
    }
}
pub fn initial(st: &StrataStack, id: u64) -> StackProvenance {
    StackProvenance {
        basement: single(
            st.basement_m,
            Setting::Initial,
            id,
            "initial-basement",
            st.basement_lith,
            st.basement_age_my,
            st.basement_age_my,
        ),
        cover: st
            .cover
            .iter()
            .enumerate()
            .map(|(i, l)| {
                single(
                    l.thickness_m,
                    Setting::Initial,
                    id,
                    &format!("initial-cover-{i}"),
                    l.lith,
                    l.age_base_my,
                    l.age_top_my,
                )
            })
            .collect(),
    }
}
pub fn validate(p: &StackProvenance, st: &StrataStack) -> Result<(), String> {
    if p.cover.len() != st.cover.len() {
        return Err("provenance/cover slot mismatch".into());
    }
    p.basement.validate(st.basement_m)?;
    for (p, l) in p.cover.iter().zip(&st.cover) {
        p.validate(l.thickness_m)?;
    }
    Ok(())
}
fn scale(p: &mut LayerProvenance, h: f64) {
    let old = p.total();
    if h == 0. {
        p.contributions.clear();
        return;
    }
    assert!(old > 0., "cannot stretch an empty provenance mixture");
    for c in &mut p.contributions {
        c.thickness_m *= h / old;
    }
}
fn combine(to: &mut LayerProvenance, from: LayerProvenance) {
    for c in from.contributions {
        let prefix = c.episode.rsplit_once(':').unwrap().0;
        if let Some(old) = to.contributions.iter_mut().find(|o| {
            o.setting == c.setting
                && o.protolith == c.protolith
                && (o.episode == c.episode
                    || (o.end_my == c.start_my && o.episode.rsplit_once(':').unwrap().0 == prefix))
        }) {
            old.thickness_m += c.thickness_m;
            old.start_my = old.start_my.min(c.start_my);
            old.end_my = old.end_my.max(c.end_my);
            old.peak_grade = old.peak_grade.max(c.peak_grade);
            if let Some(v) = c.peak_burial_m {
                old.peak_burial_m = Some(old.peak_burial_m.unwrap_or(0.).max(v));
            }
        } else {
            to.contributions.push(c);
        }
    }
}
fn slot(p: &StackProvenance, i: u32) -> &LayerProvenance {
    if i == BASEMENT {
        &p.basement
    } else {
        &p.cover[i as usize]
    }
}
fn slot_mut(p: &mut StackProvenance, i: u32) -> &mut LayerProvenance {
    if i == BASEMENT {
        &mut p.basement
    } else {
        &mut p.cover[i as usize]
    }
}
fn height(s: &StrataStack, i: u32) -> f64 {
    if i == BASEMENT {
        s.basement_m
    } else {
        s.cover[i as usize].thickness_m
    }
}

pub struct Recorder {
    pub state: StackProvenance,
    pub carrier: u64,
    pub amounts: ParcelAmounts,
    pub params: StrataParams,
    pub t0: f64,
    pub t1: f64,
}
impl Recorder {
    pub fn capture(&mut self, c: &LayerChange) {
        crate::contribution_capture::record(self, c);
        let mut next = StackProvenance {
            basement: LayerProvenance::default(),
            cover: vec![LayerProvenance::default(); c.after.cover.len()],
        };
        // Each normal link transfers the whole pre-operation mixture. The sole
        // split is the native pluton roof transfer out of basement (kind 7).
        let roof_target = c
            .links
            .iter()
            .find(|[i, j]| *i == BASEMENT && *j != BASEMENT)
            .map(|x| x[1]);
        let roof = roof_target.map_or(0., |j| {
            height(&c.after, j)
                - c.links
                    .iter()
                    .filter(|[i, k]| *i != BASEMENT && *k == j)
                    .map(|[i, _]| height(&c.before, *i))
                    .sum::<f64>()
        });
        for &[i, j] in &c.links {
            let mut part = slot(&self.state, i).clone();
            if i == BASEMENT && roof_target.is_some() {
                let h = if j == BASEMENT {
                    (c.before.basement_m - roof).max(0.)
                } else {
                    roof
                };
                if part.total() > 0. {
                    scale(&mut part, h)
                }
            }
            combine(slot_mut(&mut next, j), part);
        }
        for j in std::iter::once(BASEMENT).chain(0..c.after.cover.len() as u32) {
            let target = height(&c.after, j).max(0.);
            let p = slot_mut(&mut next, j);
            let old = p.total();
            // Only native source-addition operations may introduce contributions.
            // Merge/split sums can differ by roundoff; they must not manufacture
            // a new unknown source. Their already transported mixture is scaled.
            if target < old || !matches!(c.kind, 6 | 8 | 9 | 10 | 12 | 13 | 14) {
                if old > 0. {
                    scale(p, target);
                } else if target > 0. {
                    combine(
                        p,
                        single(
                            target,
                            Setting::Unknown,
                            self.carrier,
                            "empty-shortening",
                            LITH_BASEMENT,
                            self.t0,
                            self.t1,
                        ),
                    );
                }
            } else if target > old {
                let dh = target - old;
                let a = &self.amounts;
                let sp = &self.params;
                let parts: Vec<(f64, Setting, u8, &str)> = match c.kind {
                    8 => vec![
                        (
                            a.arc * (1. - sp.arc_extrusive_fraction),
                            Setting::Arc,
                            LITH_GRANITE,
                            "arc-pluton",
                        ),
                        (
                            a.short_pos * sp.shortening_pluton_fraction,
                            Setting::Collision,
                            LITH_GRANITE,
                            "collision-pluton",
                        ),
                    ],
                    9 => vec![(1., Setting::Plume, LITH_OCEANIC, "plume-underplate")],
                    12 => vec![
                        (
                            a.arc * sp.arc_extrusive_fraction,
                            Setting::Arc,
                            LITH_BASALT,
                            "arc-volcanic",
                        ),
                        (
                            a.plume * sp.plume_extrusive_fraction,
                            Setting::Plume,
                            LITH_BASALT,
                            "plume-volcanic",
                        ),
                    ],
                    13 => vec![(
                        1.,
                        Setting::Sedimentary,
                        DEPOSIT_CLASSES[a.class as usize].0,
                        "deposition",
                    )],
                    _ => vec![(
                        1.,
                        Setting::Unknown,
                        LITH_BASEMENT,
                        "unresolved-basement-flow",
                    )],
                };
                let sum: f64 = parts.iter().map(|p| p.0).sum();
                assert!(sum > 0., "formation addition lacks source");
                for (weight, setting, lith, process) in parts {
                    if weight > 0. {
                        combine(
                            p,
                            single(
                                dh * (weight / sum),
                                setting,
                                self.carrier,
                                process,
                                lith,
                                self.t0,
                                self.t1,
                            ),
                        );
                    }
                }
            }
            if c.kind == 18 && j != BASEMENT {
                let l = &c.after.cover[j as usize];
                for part in &mut p.contributions {
                    part.peak_grade = Some(part.peak_grade.unwrap_or(0).max(l.grade));
                    part.peak_burial_m = Some(part.peak_burial_m.unwrap_or(0.).max(l.max_burial_m));
                }
            }
            // Signed negative basement exists transiently during underflow. It is
            // a deficit, not material. The following mutation takes it from cover.
            p.validate(target)
                .expect("native formation contribution closure");
        }
        self.state = next;
    }
}
