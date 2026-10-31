//! External, immutable observation of native carriers at each step boundary.
//! No simulation fields, checkpoint codecs, or physical source hashes change.
//! Positions are rounded to microdegrees (<= 0.5e-6 degree per coordinate).
//! A removed carrier has no later position; retirement is not extrapolation.
use crate::conserve::ParcelLayer;
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs::{File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::Path,
};
type Result<T> = std::result::Result<T, String>;
const MAGIC: &[u8; 8] = b"PFTRAJ01";
pub const SCALE: f64 = 1_000_000.;
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample {
    pub step: u32,
    pub lat_microdeg: i32,
    pub lon_microdeg: i32,
    pub plate: u32,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    pub step: u32,
    pub time_my: f64,
    pub carriers: u64,
}
pub struct Capture {
    writer: BufWriter<File>,
    pub frames: Vec<Frame>,
    pub records: u64,
    pub max_error_deg: f64,
}
fn io<T>(r: std::io::Result<T>) -> Result<T> {
    r.map_err(|e| e.to_string())
}
impl Capture {
    pub fn start(path: &Path) -> Result<Self> {
        let mut writer = BufWriter::new(io(OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path))?);
        io(writer.write_all(MAGIC))?;
        Ok(Self {
            writer,
            frames: vec![],
            records: 0,
            max_error_deg: 0.,
        })
    }
    pub fn observe(&mut self, step: u32, time_my: f64, p: &ParcelLayer) -> Result<()> {
        if step as usize != self.frames.len()
            || !time_my.is_finite()
            || self.frames.last().is_some_and(|f| time_my <= f.time_my)
        {
            return Err("trajectory clock is not consecutive".into());
        }
        let ids = p
            .identity
            .as_ref()
            .ok_or("trajectory requires native identity")?
            .ids();
        if ids.len() != p.pos.len()
            || ids.len() != p.plate.len()
            || ids.windows(2).any(|w| w[0] >= w[1])
        {
            return Err("trajectory carrier shape/order".into());
        }
        io(self.writer.write_all(&step.to_le_bytes()))?;
        io(self.writer.write_all(&time_my.to_le_bytes()))?;
        io(self.writer.write_all(&(ids.len() as u64).to_le_bytes()))?;
        for (i, &id) in ids.iter().enumerate() {
            let [x, y, z] = p.pos[i];
            let norm = (x * x + y * y + z * z).sqrt();
            if !norm.is_finite() || (norm - 1.).abs() > 1e-10 {
                return Err("trajectory non-unit position".into());
            }
            let lat = z.atan2(x.hypot(y)).to_degrees();
            let lon = y.atan2(x).to_degrees();
            let qlat = (lat * SCALE).round() as i32;
            let qlon = (lon * SCALE).round() as i32;
            self.max_error_deg = self
                .max_error_deg
                .max((lat - qlat as f64 / SCALE).abs())
                .max((lon - qlon as f64 / SCALE).abs());
            for bytes in [
                &id.to_le_bytes()[..],
                &qlat.to_le_bytes()[..],
                &qlon.to_le_bytes()[..],
                &p.plate[i].0.to_le_bytes()[..],
            ] {
                io(self.writer.write_all(bytes))?;
            }
        }
        self.records += ids.len() as u64;
        self.frames.push(Frame {
            step,
            time_my,
            carriers: ids.len() as u64,
        });
        Ok(())
    }
    pub fn finish(mut self) -> Result<(Vec<Frame>, u64, f64)> {
        io(self.writer.flush())?;
        Ok((self.frames, self.records, self.max_error_deg))
    }
}
#[derive(Debug)]
pub struct Archive {
    pub frames: Vec<Frame>,
    pub tracks: BTreeMap<u64, Vec<Sample>>,
    pub records: u64,
}
impl Archive {
    pub fn read(path: &Path) -> Result<Self> {
        let len = io(std::fs::metadata(path))?.len();
        let mut r = BufReader::new(io(File::open(path))?);
        let mut magic = [0u8; 8];
        io(r.read_exact(&mut magic))?;
        if &magic != MAGIC {
            return Err("trajectory schema magic".into());
        }
        let mut a = Self {
            frames: vec![],
            tracks: BTreeMap::new(),
            records: 0,
        };
        let mut consumed = 8;
        while consumed < len {
            let mut h = [0u8; 20];
            io(r.read_exact(&mut h))?;
            consumed += 20;
            let step = u32::from_le_bytes(h[..4].try_into().unwrap());
            let time = f64::from_le_bytes(h[4..12].try_into().unwrap());
            let n = u64::from_le_bytes(h[12..20].try_into().unwrap());
            if step as usize != a.frames.len()
                || !time.is_finite()
                || a.frames.last().is_some_and(|f| time <= f.time_my)
                || n > (len - consumed) / 20
            {
                return Err("trajectory frame bounds/clock".into());
            }
            let mut prev = None;
            for _ in 0..n {
                let mut b = [0u8; 20];
                io(r.read_exact(&mut b))?;
                consumed += 20;
                let id = u64::from_le_bytes(b[..8].try_into().unwrap());
                let s = Sample {
                    step,
                    lat_microdeg: i32::from_le_bytes(b[8..12].try_into().unwrap()),
                    lon_microdeg: i32::from_le_bytes(b[12..16].try_into().unwrap()),
                    plate: u32::from_le_bytes(b[16..20].try_into().unwrap()),
                };
                if prev.is_some_and(|p| id <= p)
                    || s.lat_microdeg.unsigned_abs() > 90_000_000
                    || s.lon_microdeg.unsigned_abs() > 180_000_000
                {
                    return Err("trajectory record bounds/order".into());
                }
                prev = Some(id);
                let track = a.tracks.entry(id).or_default();
                if track.last().is_some_and(|p| p.step + 1 != step) {
                    return Err("trajectory gap/reappearing identity".into());
                }
                track.push(s);
            }
            a.records += n;
            a.frames.push(Frame {
                step,
                time_my: time,
                carriers: n,
            });
        }
        if a.frames.is_empty() {
            return Err("empty trajectory".into());
        }
        Ok(a)
    }
    pub fn sample(&self, id: u64, step: u32) -> Option<&Sample> {
        let t = self.tracks.get(&id)?;
        let first = t.first()?.step;
        t.get(step.checked_sub(first)? as usize)
            .filter(|p| p.step == step)
    }
    /// An interval must be exactly bounded by native model clock samples.
    pub fn steps(&self, start: f64, end: f64) -> Result<Vec<u32>> {
        interval_steps(&self.frames, start, end)
    }
}
pub fn interval_steps(frames: &[Frame], start: f64, end: f64) -> Result<Vec<u32>> {
    if !start.is_finite() || !end.is_finite() || start >= end {
        return Err("no simulated deposition interval".into());
    }
    let a = frames
        .iter()
        .position(|f| f.time_my == start)
        .ok_or("formation start outside native steps")?;
    let b = frames
        .iter()
        .position(|f| f.time_my == end)
        .ok_or("formation end outside native steps")?;
    if b <= a {
        return Err("reversed formation clock".into());
    }
    Ok(((a + 1) as u32..=b as u32).collect())
}
/// Every final published layer has a row. Initial/unknown and noncarrier rock
/// explicitly has no complete process-time join. Coalesced layers have many steps.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerSteps {
    pub cell: usize,
    pub layer: usize,
    pub carrier: Option<u64>,
    pub formation_steps: Vec<u32>,
    pub deposition_steps: Vec<u32>,
    pub complete: bool,
    pub reason: String,
}
pub fn layer_steps(
    cell: usize,
    layer: usize,
    carrier: Option<u64>,
    p: Option<&pf_contracts::layer_provenance::LayerProvenance>,
    frames: &[Frame],
) -> LayerSteps {
    use pf_contracts::layer_provenance::Setting;
    let mut row = LayerSteps {
        cell,
        layer,
        carrier,
        formation_steps: vec![],
        deposition_steps: vec![],
        complete: true,
        reason: "native retained contribution intervals joined to step-end positions".into(),
    };
    let Some(p) = p else {
        row.complete = false;
        row.reason = "noncarrier oceanic column has no native deposition operations".into();
        return row;
    };
    if p.contributions.is_empty() {
        row.complete = false;
        row.reason = "empty retained source".into();
    }
    for c in &p.contributions {
        if matches!(c.setting, Setting::Initial | Setting::Unknown) {
            row.complete = false;
            row.reason="initial/imported or unknown formation contribution has no recorded deposition operation".into();
            continue;
        }
        match interval_steps(frames, c.start_my, c.end_my) {
            Ok(steps) => {
                row.formation_steps.extend(&steps);
                if c.setting == Setting::Sedimentary {
                    row.deposition_steps.extend(steps);
                }
            }
            Err(e) => {
                row.complete = false;
                row.reason = e;
            }
        }
    }
    row.formation_steps.sort_unstable();
    row.formation_steps.dedup();
    row.deposition_steps.sort_unstable();
    row.deposition_steps.dedup();
    row
}

impl Frame {
    pub fn to_json(&self) -> Value {
        json!({"step":self.step,"time_my":self.time_my,"carriers":self.carriers})
    }
}
impl LayerSteps {
    pub fn to_json(&self) -> Value {
        json!({"cell":self.cell,"layer":self.layer,"carrier":self.carrier,"formation_steps":self.formation_steps,"deposition_steps":self.deposition_steps,"complete":self.complete,"reason":self.reason})
    }
    pub fn from_json(v: Value) -> Result<Self> {
        if v.as_object().is_none_or(|o| o.len() != 7) {
            return Err("formation row schema".into());
        }
        Ok(Self {
            cell: serde_json::from_value(v["cell"].clone()).map_err(|e| e.to_string())?,
            layer: serde_json::from_value(v["layer"].clone()).map_err(|e| e.to_string())?,
            carrier: serde_json::from_value(v["carrier"].clone()).map_err(|e| e.to_string())?,
            formation_steps: serde_json::from_value(v["formation_steps"].clone())
                .map_err(|e| e.to_string())?,
            deposition_steps: serde_json::from_value(v["deposition_steps"].clone())
                .map_err(|e| e.to_string())?,
            complete: v["complete"].as_bool().ok_or("formation complete")?,
            reason: v["reason"].as_str().ok_or("formation reason")?.into(),
        })
    }
}
