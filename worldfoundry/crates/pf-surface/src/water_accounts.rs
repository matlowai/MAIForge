//! Exact paired water-volume bookkeeping for native world/regional coupling.
//! One quantum is2^-32m3 (~0.233microlitres). This is arithmetic resolution,
//! not physical accuracy or proof that a routing capacity is available water.
//! Caller must validate source inventory, body connectivity and world ownership.
//! Large balances must remain integer units; f64 displays cannot retain small debits.
use pf_contracts::sim_module::{
    codec::{Reader, Writer},
    digest_bytes, SimModuleError,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
const SCALE: f64 = 4_294_967_296.;
pub const QUANTUM_M3: f64 = 1. / SCALE;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct WaterVolume(u128);
impl From<WaterVolume> for String {
    fn from(v: WaterVolume) -> Self {
        v.0.to_string()
    }
}
impl TryFrom<String> for WaterVolume {
    type Error = String;
    fn try_from(s: String) -> Result<Self, String> {
        let v = s.parse::<u128>().map_err(|_| "invalid exact water units")?;
        if v.to_string() != s {
            return Err("noncanonical exact water units".into());
        }
        Ok(Self(v))
    }
}
#[derive(Clone, Copy, Debug, Serialize)]
pub struct RoundedWaterVolume {
    pub volume: WaterVolume,
    pub requested_m3: f64,
    pub applied_m3: f64,
    /// applied minus requested
    pub rounding_m3: f64,
}
impl WaterVolume {
    pub const ZERO: Self = Self(0);
    /// Accumulate only the fractional part in f64; whole units never pass
    /// through a large running f64 balance. Caller must checkpoint remainder.
    pub fn accumulate_m3(value: f64, remainder_units: f64) -> Result<(Self, f64), String> {
        let scaled = value * SCALE;
        if !value.is_finite()
            || value < 0.
            || !scaled.is_finite()
            || scaled >= u128::MAX as f64
            || !remainder_units.is_finite()
            || !(0.0..1.0).contains(&remainder_units)
        {
            return Err("invalid fractional water accumulation".into());
        }
        let whole = scaled.floor();
        let fraction = (scaled - whole) + remainder_units;
        let carry = fraction.floor();
        let units = (whole as u128)
            .checked_add(carry as u128)
            .ok_or("water accumulation overflow")?;
        Ok((Self(units), fraction - carry))
    }
    pub fn units(self) -> u128 {
        self.0
    }
    /// Approximate display only; do not subtract world-scale f64 balances.
    pub fn approximate_m3(self) -> f64 {
        self.0 as f64 / SCALE
    }
    /// Positive sub-quantum amounts refuse rather than disappearing. Caller must
    /// retain/aggregate them. Every accepted conversion exposes its rounding.
    pub fn from_m3(v: f64) -> Result<RoundedWaterVolume, String> {
        let scaled = v * SCALE;
        if !v.is_finite() || v < 0. || !scaled.is_finite() || scaled >= u128::MAX as f64 {
            return Err("water volume outside exact-unit range".into());
        }
        let units = scaled.round() as u128;
        if v > 0. && (v < QUANTUM_M3 || units == 0) {
            return Err("positive water transfer below quantum; aggregate before admission".into());
        }
        let volume = Self(units);
        let applied = volume.approximate_m3();
        Ok(RoundedWaterVolume {
            volume,
            requested_m3: v,
            applied_m3: applied,
            rounding_m3: applied - v,
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterTransfer {
    pub id: String,
    pub source_identity: String,
    pub from: String,
    pub to: String,
    pub volume: WaterVolume,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct TransferReceipt {
    pub transfer: WaterTransfer,
    pub sequence: u64,
    pub from_before: WaterVolume,
    pub from_after: WaterVolume,
    pub to_before: WaterVolume,
    pub to_after: WaterVolume,
}
#[derive(Clone, Debug, Serialize)]
pub struct TransferOutcome {
    pub receipt: TransferReceipt,
    pub applied_now: bool,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WaterAccounts {
    source_identity: String,
    initial_account: String,
    initial_total: WaterVolume,
    accounts: BTreeMap<String, WaterVolume>,
    receipts: BTreeMap<String, TransferReceipt>,
}
fn identity(s: &str) -> Result<(), String> {
    if pf_contracts::region::is_blake3_digest(s) {
        Ok(())
    } else {
        Err("source-bound water identity required".into())
    }
}
impl WaterAccounts {
    pub fn new(
        source_identity: String,
        initial_account: String,
        volume: WaterVolume,
    ) -> Result<Self, String> {
        identity(&source_identity)?;
        identity(&initial_account)?;
        let mut accounts = BTreeMap::new();
        accounts.insert(initial_account.clone(), volume);
        Ok(Self {
            source_identity,
            initial_account,
            initial_total: volume,
            accounts,
            receipts: BTreeMap::new(),
        })
    }
    pub fn create_empty_account(&mut self, id: String) -> Result<(), String> {
        identity(&id)?;
        if self.accounts.len() >= 1_000_000 || self.accounts.contains_key(&id) {
            return Err("duplicate water account or account limit".into());
        }
        self.accounts.insert(id, WaterVolume::ZERO);
        Ok(())
    }
    pub fn balance(&self, id: &str) -> Option<WaterVolume> {
        self.accounts.get(id).copied()
    }
    pub fn source_identity(&self) -> &str {
        &self.source_identity
    }
    pub fn total(&self) -> WaterVolume {
        self.initial_total
    }
    pub fn transfer_count(&self) -> usize {
        self.receipts.len()
    }
    /// A paired transfer is atomic. An identical replay returns applied_now=false;
    /// reuse of an identity for different source/accounts/volume is an error.
    pub fn transfer(&mut self, t: WaterTransfer) -> Result<TransferOutcome, String> {
        identity(&t.id)?;
        if t.source_identity != self.source_identity {
            return Err("water transfer belongs to another source".into());
        }
        if let Some(old) = self.receipts.get(&t.id) {
            return if old.transfer == t {
                Ok(TransferOutcome {
                    receipt: old.clone(),
                    applied_now: false,
                })
            } else {
                Err("water transfer identity reused with different operation".into())
            };
        }
        if t.from == t.to || t.volume == WaterVolume::ZERO || self.receipts.len() >= 1_000_000 {
            return Err("invalid self/zero/over-limit water transfer".into());
        }
        let before_from = *self
            .accounts
            .get(&t.from)
            .ok_or("unknown water source account")?;
        let before_to = *self
            .accounts
            .get(&t.to)
            .ok_or("unknown water destination account")?;
        let after_from = WaterVolume(
            before_from
                .0
                .checked_sub(t.volume.0)
                .ok_or("water account overdraft")?,
        );
        let after_to = WaterVolume(
            before_to
                .0
                .checked_add(t.volume.0)
                .ok_or("water account overflow")?,
        );
        let receipt = TransferReceipt {
            transfer: t,
            sequence: self.receipts.len() as u64,
            from_before: before_from,
            from_after: after_from,
            to_before: before_to,
            to_after: after_to,
        };
        self.accounts
            .insert(receipt.transfer.from.clone(), after_from);
        self.accounts.insert(receipt.transfer.to.clone(), after_to);
        self.receipts
            .insert(receipt.transfer.id.clone(), receipt.clone());
        Ok(TransferOutcome {
            receipt,
            applied_now: true,
        })
    }
    /// Candidate batch: no debit survives if any later transfer refuses. Cloning is
    /// explicit O(retained accounts+receipts); single transfers are O(log n).
    pub fn transfer_batch(
        &mut self,
        transfers: &[WaterTransfer],
    ) -> Result<Vec<TransferOutcome>, String> {
        let mut candidate = self.clone();
        let mut results = Vec::with_capacity(transfers.len());
        for t in transfers {
            results.push(candidate.transfer(t.clone())?);
        }
        *self = candidate;
        Ok(results)
    }
    // Private state is created by checked constructors or fully replay-validated
    // decode_state, and all public mutations validate each new transaction before
    // commit. Encoding need not replay immutable historical receipts again.
    fn validate_balances(&self) -> Result<(), String> {
        identity(&self.source_identity)?;
        identity(&self.initial_account)?;
        if self.accounts.is_empty()
            || self.accounts.len() > 1_000_000
            || self.receipts.len() > 1_000_000
            || !self.accounts.contains_key(&self.initial_account)
        {
            return Err("invalid water ledger dimensions".into());
        }
        let mut sum = 0u128;
        for (k, v) in &self.accounts {
            identity(k)?;
            sum = sum.checked_add(v.0).ok_or("water total overflow")?;
        }
        if sum != self.initial_total.0 {
            return Err("water inventory total changed".into());
        }
        Ok(())
    }
    // Full replay is mandatory at the untrusted checkpoint boundary, even if an
    // altered payload has been rehashed. Do not replace this in decode_state.
    fn validate(&self) -> Result<(), String> {
        self.validate_balances()?;
        let mut replay = BTreeMap::new();
        for k in self.accounts.keys() {
            replay.insert(k.clone(), WaterVolume::ZERO);
        }
        replay.insert(self.initial_account.clone(), self.initial_total);
        let mut receipts: Vec<_> = self.receipts.iter().collect();
        receipts.sort_by_key(|(_, v)| v.sequence);
        for (seq, (key, r)) in receipts.into_iter().enumerate() {
            let t = &r.transfer;
            identity(&t.id)?;
            if key != &t.id
                || r.sequence != seq as u64
                || t.source_identity != self.source_identity
                || t.from == t.to
                || t.volume == WaterVolume::ZERO
            {
                return Err("water receipt identity/order mismatch".into());
            }
            let a = *replay.get(&t.from).ok_or("receipt source absent")?;
            let b = *replay.get(&t.to).ok_or("receipt destination absent")?;
            let aa = WaterVolume(a.0.checked_sub(t.volume.0).ok_or("receipt overdraft")?);
            let bb = WaterVolume(b.0.checked_add(t.volume.0).ok_or("receipt overflow")?);
            if (a, b, aa, bb) != (r.from_before, r.to_before, r.from_after, r.to_after) {
                return Err("water receipt balance mismatch".into());
            }
            replay.insert(t.from.clone(), aa);
            replay.insert(t.to.clone(), bb);
        }
        if replay != self.accounts {
            return Err("water ledger differs from receipt replay".into());
        }
        Ok(())
    }
    pub fn encode_state(&self) -> Result<Vec<u8>, String> {
        self.validate_balances()?;
        let mut w = Writer::new();
        w.tag(b"WAC1");
        w.str(&self.source_identity);
        w.str(&self.initial_account);
        put_volume(&mut w, self.initial_total);
        w.len(self.accounts.len());
        for (k, v) in &self.accounts {
            w.str(k);
            put_volume(&mut w, *v);
        }
        w.len(self.receipts.len());
        for r in self.receipts.values() {
            let t = &r.transfer;
            w.str(&t.id);
            w.str(&t.source_identity);
            w.str(&t.from);
            w.str(&t.to);
            put_volume(&mut w, t.volume);
            w.u64(r.sequence);
            for v in [r.from_before, r.from_after, r.to_before, r.to_after] {
                put_volume(&mut w, v);
            }
        }
        let mut frame = Writer::new();
        frame.tag(b"WAF1");
        frame.str(&digest_bytes(&w.buf));
        frame.bytes(&w.buf);
        Ok(frame.buf)
    }
    pub fn decode_state(bytes: &[u8], expected_source: &str) -> Result<Self, String> {
        if bytes.len() > 512 * 1024 * 1024 {
            return Err("water ledger checkpoint size limit".into());
        }
        let parse = || -> Result<Self, SimModuleError> {
            let mut f = Reader::new(bytes);
            f.tag(b"WAF1")?;
            let digest = f.str()?;
            let payload = f.bytes()?;
            f.finish()?;
            if digest_bytes(payload) != digest {
                return Err(SimModuleError::Checkpoint(
                    "water account digest mismatch".into(),
                ));
            }
            let mut r = Reader::new(payload);
            r.tag(b"WAC1")?;
            let source_identity = r.str()?;
            let initial_account = r.str()?;
            let initial_total = get_volume(&mut r)?;
            let n = r.len()?;
            if n > 1_000_000 {
                return Err(SimModuleError::Checkpoint("water account count".into()));
            }
            let mut accounts = BTreeMap::new();
            for _ in 0..n {
                let k = r.str()?;
                let v = get_volume(&mut r)?;
                if accounts.insert(k, v).is_some() {
                    return Err(SimModuleError::Checkpoint("duplicate water account".into()));
                }
            }
            let n = r.len()?;
            if n > 1_000_000 {
                return Err(SimModuleError::Checkpoint("water receipt count".into()));
            }
            let mut receipts = BTreeMap::new();
            for _ in 0..n {
                let id = r.str()?;
                let t = WaterTransfer {
                    id: id.clone(),
                    source_identity: r.str()?,
                    from: r.str()?,
                    to: r.str()?,
                    volume: get_volume(&mut r)?,
                };
                let receipt = TransferReceipt {
                    transfer: t,
                    sequence: r.u64()?,
                    from_before: get_volume(&mut r)?,
                    from_after: get_volume(&mut r)?,
                    to_before: get_volume(&mut r)?,
                    to_after: get_volume(&mut r)?,
                };
                if receipts.insert(id, receipt).is_some() {
                    return Err(SimModuleError::Checkpoint("duplicate water receipt".into()));
                }
            }
            r.finish()?;
            Ok(Self {
                source_identity,
                initial_account,
                initial_total,
                accounts,
                receipts,
            })
        };
        let state = parse().map_err(|e| e.to_string())?;
        if state.source_identity != expected_source {
            return Err("water ledger source mismatch".into());
        }
        state.validate()?;
        Ok(state)
    }
}
fn put_volume(w: &mut Writer, v: WaterVolume) {
    w.u64((v.0 >> 64) as u64);
    w.u64(v.0 as u64);
}
fn get_volume(r: &mut Reader<'_>) -> Result<WaterVolume, SimModuleError> {
    Ok(WaterVolume(((r.u64()? as u128) << 64) | r.u64()? as u128))
}

/// Convert the existing hydraulic kernel's identified boundary flows into paired
/// native transactions. This does not establish/admit either account's inventory.
/// Anonymous volumetric forcing must first obtain its own identified source account.
#[derive(Clone, Debug, Serialize)]
pub struct BoundaryTransferBatch {
    pub transfers: Vec<WaterTransfer>,
    pub incoming_rounding_m3: f64,
    pub outgoing_rounding_m3: f64,
}
pub fn boundary_transfers(
    source_identity: &str,
    step_identity: &str,
    regional_account: &str,
    report: &crate::river_water::WaterStep,
) -> Result<BoundaryTransferBatch, String> {
    identity(source_identity)?;
    identity(step_identity)?;
    identity(regional_account)?;
    if report.supplied_m3 != 0. {
        return Err("volumetric water sources require identified funding before settlement".into());
    }
    let mut seen = std::collections::BTreeSet::new();
    let mut incoming = Vec::new();
    let mut outgoing = Vec::new();
    let (mut qi, mut qo, mut ri, mut ro) = (0., 0., 0., 0.);
    for v in &report.boundary_volumes {
        identity(&v.reservoir_identity)?;
        if v.reservoir_identity == regional_account
            || !seen.insert(&v.reservoir_identity)
            || !v.incoming_m3.is_finite()
            || !v.outgoing_m3.is_finite()
            || v.incoming_m3 < 0.
            || v.outgoing_m3 < 0.
        {
            return Err("invalid native boundary transfer record".into());
        }
        qi += v.incoming_m3;
        qo += v.outgoing_m3;
        for (direction, amount) in [("in", v.incoming_m3), ("out", v.outgoing_m3)] {
            if amount == 0. {
                continue;
            }
            let rounded = WaterVolume::from_m3(amount)?;
            let id = digest_bytes(
                format!(
                    "{source_identity}:{step_identity}:{regional_account}:{}:{direction}",
                    v.reservoir_identity
                )
                .as_bytes(),
            );
            let (from, to) = if direction == "in" {
                (v.reservoir_identity.clone(), regional_account.to_owned())
            } else {
                (regional_account.to_owned(), v.reservoir_identity.clone())
            };
            let t = WaterTransfer {
                id,
                source_identity: source_identity.into(),
                from,
                to,
                volume: rounded.volume,
            };
            if direction == "in" {
                ri += rounded.rounding_m3;
                incoming.push(t);
            } else {
                ro += rounded.rounding_m3;
                outgoing.push(t);
            }
        }
    }
    if !qi.is_finite()
        || !qo.is_finite()
        || qi.to_bits() != report.boundary_in_m3.to_bits()
        || qo.to_bits() != report.boundary_out_m3.to_bits()
    {
        return Err("boundary totals do not match native hydraulic report".into());
    }
    incoming.extend(outgoing);
    Ok(BoundaryTransferBatch {
        transfers: incoming,
        incoming_rounding_m3: ri,
        outgoing_rounding_m3: ro,
    })
}
