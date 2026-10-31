//! Bounded native 2D hydraulic comparison on fixed terrain, not a meander model.
//! First-order Cartesian finite volumes, Rusanov flux, Audusse et al. (2004)
//! hydrostatic reconstruction, reflective outer walls and explicit volume sources.
//! https://publications.imp.fu-berlin.de/478/1/file_2004_siam.pdf
//! Depth-averaged hydrostatic pressure; no vertical jets, sediment or bank motion.
//! Delestre et al. warn of slope/depth/mesh limitations (arxiv.org/abs/1206.4986).
//! Moving wet faces truncated by a bed step are counted, never silently qualified.
//! Fixed-wall domains are comparison controls, not automatically world boundaries.
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct WaterCell {
    pub depth_m: f64,
    pub flux_x_m2_s: f64,
    pub flux_y_m2_s: f64,
}
impl WaterCell {
    fn velocity(self) -> (f64, f64) {
        if self.depth_m == 0. {
            (0., 0.)
        } else {
            (
                self.flux_x_m2_s / self.depth_m,
                self.flux_y_m2_s / self.depth_m,
            )
        }
    }
    fn valid(self) -> bool {
        self.depth_m.is_finite()
            && self.depth_m >= 0.
            && self.flux_x_m2_s.is_finite()
            && self.flux_y_m2_s.is_finite()
            && (self.depth_m > 0. || (self.flux_x_m2_s == 0. && self.flux_y_m2_s == 0.))
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WaterGrid {
    pub source_identity: String,
    pub nx: usize,
    pub ny: usize,
    pub spacing_m: f64,
    pub gravity_m_s2: f64,
    pub friction_cf: f64,
    pub bed_m: Vec<f64>,
    pub water: Vec<WaterCell>,
    pub seconds: f64,
    pub steps: u64,
}
/// Resolved external volume and advected momentum; not an inferred river inlet.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WaterSource {
    pub cell: usize,
    pub volume_m3_s: f64,
    pub velocity_x_m_s: f64,
    pub velocity_y_m_s: f64,
}
/// Prescribed exterior reservoir state at an outer face. Not a free radiation
/// condition or proof of the correct stage. Every exchanged volume is reported.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Edge {
    West,
    East,
    South,
    North,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WaterBoundary {
    pub edge: Edge,
    pub index: usize,
    pub reservoir_identity: String,
    pub exterior_bed_m: f64,
    pub exterior_water: WaterCell,
}
#[derive(Clone, Debug, Serialize)]
pub struct BoundaryVolume {
    pub reservoir_identity: String,
    pub incoming_m3: f64,
    pub outgoing_m3: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct WaterStep {
    pub dt_seconds: f64,
    pub initial_m3: f64,
    pub supplied_m3: f64,
    pub boundary_in_m3: f64,
    pub boundary_out_m3: f64,
    pub boundary_volumes: Vec<BoundaryVolume>,
    pub final_m3: f64,
    pub residual_m3: f64,
    pub moving_wet_faces_truncated_by_bed: usize,
    pub max_froude: f64,
    pub max_speed_m_s: f64,
}
impl WaterGrid {
    pub fn validate(&self) -> Result<(), String> {
        let n = self.nx.checked_mul(self.ny).ok_or("grid size overflow")?;
        if !pf_contracts::region::is_blake3_digest(&self.source_identity)
            || self.nx < 2
            || self.ny < 2
            || n > 1_048_576
            || self.bed_m.len() != n
            || self.water.len() != n
            || !self.spacing_m.is_finite()
            || self.spacing_m <= 0.
            || !(self.spacing_m * self.spacing_m).is_finite()
            || self.spacing_m * self.spacing_m == 0.
            || !self.gravity_m_s2.is_finite()
            || self.gravity_m_s2 <= 0.
            || !self.friction_cf.is_finite()
            || self.friction_cf < 0.
            || !self.seconds.is_finite()
            || self.seconds < 0.
            || self.bed_m.iter().any(|z| !z.is_finite())
            || self.water.iter().any(|c| !c.valid())
        {
            return Err("invalid bounded river-water state".into());
        }
        if !self.volume_m3().is_finite() {
            return Err("water volume overflow".into());
        }
        Ok(())
    }
    /// Bit-exact native state payload with corruption detection. The host must also
    /// bind physical forcing/schedule to its recipe; this is not world admission.
    pub fn encode_state(&self) -> Result<Vec<u8>, String> {
        use pf_contracts::sim_module::{codec::Writer, digest_bytes};
        self.validate()?;
        let mut p = Writer::new();
        p.str(&self.source_identity);
        p.u64(self.nx as u64);
        p.u64(self.ny as u64);
        for x in [
            self.spacing_m,
            self.gravity_m_s2,
            self.friction_cf,
            self.seconds,
        ] {
            p.f64(x);
        }
        p.u64(self.steps);
        p.f64s(&self.bed_m);
        for c in &self.water {
            p.f64(c.depth_m);
            p.f64(c.flux_x_m2_s);
            p.f64(c.flux_y_m2_s);
        }
        let mut w = Writer::new();
        w.tag(b"RWC1");
        w.str(&digest_bytes(&p.buf));
        w.bytes(&p.buf);
        Ok(w.buf)
    }
    pub fn decode_state(bytes: &[u8], expected_source: &str) -> Result<Self, String> {
        use pf_contracts::sim_module::{codec::Reader, digest_bytes, SimModuleError};
        if bytes.len() > 40 * 1024 * 1024 {
            return Err("water checkpoint exceeds bounded grid".into());
        }
        let parse = || -> Result<Self, SimModuleError> {
            let mut frame = Reader::new(bytes);
            frame.tag(b"RWC1")?;
            let digest = frame.str()?;
            let payload = frame.bytes()?;
            frame.finish()?;
            if digest_bytes(payload) != digest {
                return Err(SimModuleError::Checkpoint(
                    "water payload digest mismatch".into(),
                ));
            }
            let mut p = Reader::new(payload);
            let source_identity = p.str()?;
            if source_identity != expected_source {
                return Err(SimModuleError::Checkpoint("water source mismatch".into()));
            }
            let nx = p.u64()?;
            let ny = p.u64()?;
            let n = nx
                .checked_mul(ny)
                .filter(|n| *n <= 1_048_576)
                .ok_or_else(|| SimModuleError::Checkpoint("water grid bound".into()))?
                as usize;
            let spacing_m = p.f64()?;
            let gravity_m_s2 = p.f64()?;
            let friction_cf = p.f64()?;
            let seconds = p.f64()?;
            let steps = p.u64()?;
            let bed_m = p.f64s()?;
            if bed_m.len() != n {
                return Err(SimModuleError::Checkpoint("water bed length".into()));
            }
            let mut water = Vec::with_capacity(n);
            for _ in 0..n {
                water.push(WaterCell {
                    depth_m: p.f64()?,
                    flux_x_m2_s: p.f64()?,
                    flux_y_m2_s: p.f64()?,
                });
            }
            p.finish()?;
            Ok(Self {
                source_identity,
                nx: nx as usize,
                ny: ny as usize,
                spacing_m,
                gravity_m_s2,
                friction_cf,
                bed_m,
                water,
                seconds,
                steps,
            })
        };
        let state = parse().map_err(|e| e.to_string())?;
        state.validate()?;
        Ok(state)
    }
    pub fn volume_m3(&self) -> f64 {
        self.water.iter().map(|c| c.depth_m).sum::<f64>() * self.spacing_m * self.spacing_m
    }
    /// Conservative unsplit 2D CFL bound. A source may require a smaller step.
    pub fn timestep_limit(&self) -> Result<f64, String> {
        self.validate()?;
        limit(self.water.iter(), self.spacing_m, self.gravity_m_s2)
    }
    fn boundary_lookup<'a>(
        &self,
        boundaries: &'a [WaterBoundary],
    ) -> Result<Vec<Option<&'a WaterBoundary>>, String> {
        let mut lookup = vec![None; 2 * (self.nx + self.ny)];
        for b in boundaries {
            let (length, offset) = match b.edge {
                Edge::West => (self.ny, 0),
                Edge::East => (self.ny, self.ny),
                Edge::South => (self.nx, 2 * self.ny),
                Edge::North => (self.nx, 2 * self.ny + self.nx),
            };
            if b.index >= length
                || !b.exterior_bed_m.is_finite()
                || !b.exterior_water.valid()
                || !pf_contracts::region::is_blake3_digest(&b.reservoir_identity)
            {
                return Err("invalid explicit water boundary".into());
            }
            if lookup[offset + b.index].replace(b).is_some() {
                return Err("duplicate water boundary face".into());
            }
        }
        Ok(lookup)
    }
    pub fn timestep_limit_open(&self, boundaries: &[WaterBoundary]) -> Result<f64, String> {
        self.validate()?;
        self.boundary_lookup(boundaries)?;
        limit(
            self.water
                .iter()
                .chain(boundaries.iter().map(|b| &b.exterior_water)),
            self.spacing_m,
            self.gravity_m_s2,
        )
    }
    /// Atomic candidate: failure never mutates water, terrain or clock. No depth clipping.
    pub fn step(&mut self, dt: f64, sources: &[WaterSource]) -> Result<WaterStep, String> {
        self.step_open(dt, sources, &[])
    }
    pub fn step_open(
        &mut self,
        dt: f64,
        sources: &[WaterSource],
        boundaries: &[WaterBoundary],
    ) -> Result<WaterStep, String> {
        self.validate()?;
        let lookup = self.boundary_lookup(boundaries)?;
        let mut transfers: std::collections::BTreeMap<String, (f64, f64)> =
            std::collections::BTreeMap::new();
        if !dt.is_finite()
            || dt <= 0.
            || !(self.seconds + dt).is_finite()
            || self.seconds + dt <= self.seconds
            || self.steps == u64::MAX
        {
            return Err("invalid river-water timestep/clock".into());
        }
        let initial = self.volume_m3();
        let area = self.spacing_m * self.spacing_m;
        let mut input = self.water.clone();
        let mut supplied = 0.;
        for s in sources {
            if s.cell >= input.len()
                || !s.volume_m3_s.is_finite()
                || s.volume_m3_s < 0.
                || !s.velocity_x_m_s.is_finite()
                || !s.velocity_y_m_s.is_finite()
            {
                return Err("invalid explicit water source".into());
            }
            let v = s.volume_m3_s * dt;
            supplied += v;
            let dh = v / area;
            input[s.cell].depth_m += dh;
            input[s.cell].flux_x_m2_s += dh * s.velocity_x_m_s;
            input[s.cell].flux_y_m2_s += dh * s.velocity_y_m_s;
        }
        if !supplied.is_finite() || input.iter().any(|c| !c.valid()) {
            return Err("water source overflow".into());
        }
        if dt
            > limit(
                input
                    .iter()
                    .chain(boundaries.iter().map(|b| &b.exterior_water)),
                self.spacing_m,
                self.gravity_m_s2,
            )?
        {
            return Err("river-water CFL timestep exceeded (including sources)".into());
        }
        let mut change = vec![[0.; 3]; input.len()];
        let factor = dt / self.spacing_m;
        let mut truncated = 0;
        for axis in 0..2 {
            let (normal, tangent) = if axis == 0 {
                (self.nx, self.ny)
            } else {
                (self.ny, self.nx)
            };
            let idx = |a: usize, b: usize| {
                if axis == 0 {
                    b * self.nx + a
                } else {
                    a * self.nx + b
                }
            };
            for b in 0..tangent {
                for face in 0..=normal {
                    let li = if face > 0 {
                        Some(idx(face - 1, b))
                    } else {
                        None
                    };
                    let ri = if face < normal {
                        Some(idx(face, b))
                    } else {
                        None
                    };
                    let reflect = |mut c: WaterCell| {
                        if axis == 0 {
                            c.flux_x_m2_s = -c.flux_x_m2_s;
                        } else {
                            c.flux_y_m2_s = -c.flux_y_m2_s;
                        }
                        c
                    };
                    let bc = if face == 0 {
                        lookup[if axis == 0 { b } else { 2 * self.ny + b }]
                    } else if face == normal {
                        lookup[if axis == 0 {
                            self.ny + b
                        } else {
                            2 * self.ny + self.nx + b
                        }]
                    } else {
                        None
                    };
                    let l = li.map(|i| input[i]).unwrap_or_else(|| {
                        bc.map(|c| c.exterior_water)
                            .unwrap_or_else(|| reflect(input[ri.unwrap()]))
                    });
                    let rr = ri.map(|i| input[i]).unwrap_or_else(|| {
                        bc.map(|c| c.exterior_water)
                            .unwrap_or_else(|| reflect(input[li.unwrap()]))
                    });
                    let zl = li.map(|i| self.bed_m[i]).unwrap_or_else(|| {
                        bc.map(|c| c.exterior_bed_m)
                            .unwrap_or(self.bed_m[ri.unwrap()])
                    });
                    let zr = ri.map(|i| self.bed_m[i]).unwrap_or_else(|| {
                        bc.map(|c| c.exterior_bed_m)
                            .unwrap_or(self.bed_m[li.unwrap()])
                    });
                    let (fl, fr, cut) = flux(l, rr, zl, zr, axis, self.gravity_m_s2);
                    if let Some(bc) = bc {
                        let incoming =
                            dt * self.spacing_m * if li.is_none() { fr[0] } else { -fl[0] };
                        if !incoming.is_finite() {
                            return Err("boundary flux overflow".into());
                        }
                        let entry = transfers.entry(bc.reservoir_identity.clone()).or_default();
                        if incoming >= 0. {
                            entry.0 += incoming;
                        } else {
                            entry.1 -= incoming;
                        }
                    }
                    if li.is_some() && ri.is_some() && cut {
                        truncated += 1;
                    }
                    for k in 0..3 {
                        if let Some(i) = li {
                            change[i][k] -= factor * fl[k];
                        }
                        if let Some(i) = ri {
                            change[i][k] += factor * fr[k];
                        }
                    }
                }
            }
        }
        let mut max_fr = 0f64;
        let mut max_speed = 0f64;
        for (c, d) in input.iter_mut().zip(change) {
            c.depth_m += d[0];
            c.flux_x_m2_s += d[1];
            c.flux_y_m2_s += d[2];
            if !c.valid() {
                return Err(
                    "nonpositive/nonfinite hydraulic candidate; no clipping performed".into(),
                );
            }
            let (u, v) = c.velocity();
            let speed = u.hypot(v);
            if c.depth_m > 0. {
                // Backward drag magnitude with explicit direction: no mass change.
                let damping = 1. + dt * self.friction_cf * speed / c.depth_m;
                if !damping.is_finite() {
                    return Err("friction arithmetic overflow".into());
                }
                c.flux_x_m2_s /= damping;
                c.flux_y_m2_s /= damping;
                let speed = speed / damping;
                max_speed = max_speed.max(speed);
                max_fr = max_fr.max(speed / (self.gravity_m_s2 * c.depth_m).sqrt());
            }
        }
        let final_m3 = input.iter().map(|c| c.depth_m).sum::<f64>() * area;
        let boundary_in_m3 = transfers.values().map(|v| v.0).sum::<f64>();
        let boundary_out_m3 = transfers.values().map(|v| v.1).sum::<f64>();
        let residual = initial + supplied + boundary_in_m3 - boundary_out_m3 - final_m3;
        if !residual.is_finite()
            || !boundary_in_m3.is_finite()
            || !boundary_out_m3.is_finite()
            || !final_m3.is_finite()
            || !max_fr.is_finite()
            || residual.abs()
                > 1e-8 + 1e-10 * (initial + supplied + boundary_in_m3 + boundary_out_m3 + final_m3)
        {
            return Err("native 2D water volume/arithmetic failed".into());
        }
        self.water = input;
        self.seconds += dt;
        self.steps += 1;
        Ok(WaterStep {
            dt_seconds: dt,
            initial_m3: initial,
            supplied_m3: supplied,
            boundary_in_m3,
            boundary_out_m3,
            boundary_volumes: transfers
                .into_iter()
                .map(
                    |(reservoir_identity, (incoming_m3, outgoing_m3))| BoundaryVolume {
                        reservoir_identity,
                        incoming_m3,
                        outgoing_m3,
                    },
                )
                .collect(),
            final_m3,
            residual_m3: residual,
            moving_wet_faces_truncated_by_bed: truncated,
            max_froude: max_fr,
            max_speed_m_s: max_speed,
        })
    }
}
fn limit<'a>(cells: impl Iterator<Item = &'a WaterCell>, dx: f64, g: f64) -> Result<f64, String> {
    let mut ax = 0f64;
    let mut ay = 0f64;
    for c in cells {
        let (u, v) = c.velocity();
        let wave = (g * c.depth_m).sqrt();
        ax = ax.max(u.abs() + wave);
        ay = ay.max(v.abs() + wave);
    }
    let speed = ax + ay;
    if !speed.is_finite() {
        return Err("wave speed overflow".into());
    }
    // Finite dry-state cap only limits execution, not scientific history.
    Ok(if speed == 0. {
        f64::MAX
    } else {
        0.2 * dx / speed
    })
}
fn flux(
    l: WaterCell,
    r: WaterCell,
    zl: f64,
    zr: f64,
    axis: usize,
    g: f64,
) -> ([f64; 3], [f64; 3], bool) {
    let (ul, vl) = l.velocity();
    let (ur, vr) = r.velocity();
    let top = zl.max(zr);
    let hl = (l.depth_m + zl - top).max(0.);
    let hr = (r.depth_m + zr - top).max(0.);
    let left = [hl, hl * ul, hl * vl];
    let right = [hr, hr * ur, hr * vr];
    let un_l = if axis == 0 { ul } else { vl };
    let un_r = if axis == 0 { ur } else { vr };
    let a = (un_l.abs() + (g * hl).sqrt()).max(un_r.abs() + (g * hr).sqrt());
    let mut f = [0.; 3];
    for k in 0..3 {
        let pl = if k == axis + 1 { 0.5 * g * hl * hl } else { 0. };
        let pr = if k == axis + 1 { 0.5 * g * hr * hr } else { 0. };
        f[k] = 0.5 * (left[k] * un_l + pl + right[k] * un_r + pr) - 0.5 * a * (right[k] - left[k]);
    }
    let mut fl = f;
    let mut fr = f;
    fl[axis + 1] += 0.5 * g * (l.depth_m * l.depth_m - hl * hl);
    fr[axis + 1] += 0.5 * g * (r.depth_m * r.depth_m - hr * hr);
    let moving = ul.abs() + vl.abs() + ur.abs() + vr.abs() > 1e-10;
    let cut = moving && l.depth_m > 0. && r.depth_m > 0. && (hl == 0. || hr == 0.);
    (fl, fr, cut)
}
