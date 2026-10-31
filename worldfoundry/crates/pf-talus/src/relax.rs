//! Repose relaxation of the deposited debris (Jacobi, synchronous, mass conserving).
//!
//! Each iteration: every post i with debris computes the excess drop to each of its 8 neighbours,
//! e_ij = (h_i - h_j) - tan(repose) d_ij (h = routing surface + debris), and releases
//! q_ij = f_i e_ij with f_i = beta / (n_i + 1) (n_i lower neighbours; equalizes every lower pair at once for equal
//! areas), scaled down so the total never exceeds T_i. A receiver recomputes
//! the same e_ij from the same bytes, so the flux leaving i equals the flux arriving at j (to rounding) and the
//! update of a post depends only on posts within 2 of it: N iterations have an influence radius of 2N posts, which
//! is what the halo budgets. Water-sink posts absorb what reaches them (counted as lost). Only debris moves; the rock
//! surface never changes. Tracers (moment fields) move with the debris in proportion.

use crate::chute::NB;

pub struct RelaxOut {
    pub lost_m3: f64,
    /// Gross volume moved between posts (sum of every release; coarse levels at their cell area, plus the volume
    /// the prolongation redistributes inside cells).
    pub moved_m3: f64,
    pub iters: u32,
}

#[allow(clippy::too_many_arguments)]
pub fn relax(nx: usize, ny: usize, dx: f64, base: &[f64], t: &mut Vec<f64>, tracers: &mut [Vec<f64>], sink: &[u8], tan_r: f64, iters: u32, beta: f64, threads: usize) -> RelaxOut {
    relax_area(nx, ny, dx, base, None, None, t, tracers, sink, tan_r, iters, beta, threads)
}

/// `relax` with a relative plan area per post (coarse cells cut by the domain edge hold fewer fine posts): a
/// release of thickness q from post i is volume q a_i, and arrives at j as thickness q a_i / a_j, so volume is
/// conserved for unequal cells.
#[allow(clippy::too_many_arguments)]
/// `tan_post` (optional): a per-post repose tangent; the SENDER's value sets the excess of a pair, so sender and
/// receiver compute the same bytes.
pub fn relax_area(nx: usize, ny: usize, dx: f64, base: &[f64], area: Option<&[f64]>, tan_post: Option<&[f64]>, t: &mut Vec<f64>, tracers: &mut [Vec<f64>], sink: &[u8], tan_r: f64, iters: u32, beta: f64, threads: usize) -> RelaxOut {
    let ones;
    let area: &[f64] = match area { Some(a) => a, None => { ones = vec![1.0f64; nx * ny]; &ones } };
    let tconst;
    let tanp: &[f64] = match tan_post { Some(a) => a, None => { tconst = vec![tan_r; nx * ny]; &tconst } };
    let n = nx * ny;
    let mut q = vec![0.0f64; n];
    let mut e_sum = vec![0.0f64; n];
    let mut t2 = vec![0.0f64; n];
    let mut tr2: Vec<Vec<f64>> = tracers.iter().map(|_| vec![0.0f64; n]).collect();
    let mut lost = 0.0f64;
    let dd: [f64; 8] = core::array::from_fn(|m| if NB[m].0 != 0 && NB[m].1 != 0 { dx * std::f64::consts::SQRT_2 } else { dx });
    let threads = threads.max(1);
    let rows_per = ny.div_ceil(threads);
    let mut lost_rows = vec![0.0f64; ny];
    let mut moved = 0.0f64;
    for _ in 0..iters {
        // pass 1: release per post
        {
            let tt: &Vec<f64> = t;
            std::thread::scope(|sc| {
                for (ci, (qc, ec)) in q.chunks_mut(rows_per * nx).zip(e_sum.chunks_mut(rows_per * nx)).enumerate() {
                    let j0 = ci * rows_per;
                    sc.spawn(move || {
                        for (r, (qr, er)) in qc.chunks_mut(nx).zip(ec.chunks_mut(nx)).enumerate() {
                            let j = j0 + r;
                            for i in 0..nx {
                                let k = j * nx + i;
                                qr[i] = 0.0;
                                er[i] = 0.0;
                                if tt[k] <= 0.0 || sink[k] != 0 { continue; }
                                let hi = base[k] + tt[k];
                                let (mut es, mut nl) = (0.0f64, 0u32);
                                for (m, &(di, dj)) in NB.iter().enumerate() {
                                    let (ni, nj) = (i as i64 + di, j as i64 + dj);
                                    if ni < 0 || nj < 0 || ni >= nx as i64 || nj >= ny as i64 { continue; }
                                    let kn = nj as usize * nx + ni as usize;
                                    let e = (hi - (base[kn] + tt[kn])) - tanp[k] * dd[m];
                                    if e > 0.0 { es += e; nl += 1; }
                                }
                                if es > 0.0 {
                                    // per-receiver share e_ij / (n + 1) equalizes every lower pair at once (equal areas);
                                    // beta damps simultaneous donors; capped by the debris present
                                    let f0 = beta / (nl as f64 + 1.0);
                                    let f = if f0 * es > tt[k] { tt[k] / es } else { f0 };
                                    qr[i] = f * es;
                                    er[i] = f;
                                }
                            }
                        }
                    });
                }
            });
        }
        moved += q.iter().zip(area.iter()).map(|(a, b)| a * b).sum::<f64>() * dx * dx;
        // pass 2: update
        {
            let tt: &Vec<f64> = t;
            let trs: &[Vec<f64>] = tracers;
            let (qq, ee) = (&q, &e_sum);
            let mut tr_chunks: Vec<Vec<&mut [f64]>> = (0..ny.div_ceil(rows_per)).map(|_| Vec::new()).collect();
            for trv in tr2.iter_mut() {
                for (ci, c) in trv.chunks_mut(rows_per * nx).enumerate() { tr_chunks[ci].push(c); }
            }
            std::thread::scope(|sc| {
                for (ci, ((tc, mut trc), lc)) in t2.chunks_mut(rows_per * nx).zip(tr_chunks.into_iter()).zip(lost_rows.chunks_mut(rows_per)).enumerate() {
                    let j0 = ci * rows_per;
                    sc.spawn(move || {
                        let rows = tc.len() / nx;
                        for r in 0..rows {
                            let j = j0 + r;
                            lc[r] = 0.0;
                            for i in 0..nx {
                                let k = j * nx + i;
                                let hk = base[k] + tt[k];
                                let mut v = tt[k] - qq[k];
                                let keep = if tt[k] > 0.0 { 1.0 - qq[k] / tt[k] } else { 0.0 };
                                let mut trv: [f64; 4] = [0.0; 4];
                                for (x, tr) in trs.iter().enumerate() { trv[x] = tr[k] * keep; }
                                for (m, &(di, dj)) in NB.iter().enumerate() {
                                    let (ni, nj) = (i as i64 + di, j as i64 + dj);
                                    if ni < 0 || nj < 0 || ni >= nx as i64 || nj >= ny as i64 { continue; }
                                    let kn = nj as usize * nx + ni as usize;
                                    if qq[kn] <= 0.0 { continue; }
                                    let e = ((base[kn] + tt[kn]) - hk) - tanp[kn] * dd[m];
                                    if e > 0.0 {
                                        let fl = ee[kn] * e;
                                        let ar = area[kn] / area[k];
                                        v += fl * ar;
                                        let share = fl / tt[kn] * ar;
                                        for (x, tr) in trs.iter().enumerate() { trv[x] += tr[kn] * share; }
                                    }
                                }
                                if sink[k] != 0 {
                                    lc[r] += v.max(0.0) * area[k];
                                    v = 0.0;
                                    trv = [0.0; 4];
                                }
                                tc[r * nx + i] = v.max(0.0);
                                for (x, c) in trc.iter_mut().enumerate() { c[r * nx + i] = trv[x]; }
                            }
                        }
                    });
                }
            });
        }
        std::mem::swap(t, &mut t2);
        for (a, b) in tracers.iter_mut().zip(tr2.iter_mut()) { std::mem::swap(a, b); }
        for l in &lost_rows { lost += l * dx * dx; }
    }
    RelaxOut { lost_m3: lost, moved_m3: moved, iters }
}

/// Share of debris posts (T > t_min) whose steepest descent on the debris surface exceeds tan(repose + tol).
pub fn residual(nx: usize, ny: usize, dx: f64, base: &[f64], t: &[f64], tan_lim: f64, t_min: f64) -> (u64, u64) {
    let (mut over, mut all) = (0u64, 0u64);
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if t[k] <= t_min { continue; }
            all += 1;
            let hk = base[k] + t[k];
            let mut smax = 0.0f64;
            for &(di, dj) in NB.iter() {
                let (ni, nj) = (i as i64 + di, j as i64 + dj);
                if ni < 0 || nj < 0 || ni >= nx as i64 || nj >= ny as i64 { continue; }
                let kn = nj as usize * nx + ni as usize;
                let d = if di != 0 && dj != 0 { dx * std::f64::consts::SQRT_2 } else { dx };
                let s = (hk - (base[kn] + t[kn])) / d;
                if s > smax { smax = s; }
            }
            if smax > tan_lim { over += 1; }
        }
    }
    (over, all)
}

/// Multilevel acceleration (coarse to fine). Level factor F aggregates F x F posts aligned on the GLOBAL lattice
/// (cell = floor(g / F)), so the same cells exist in every tile domain (cells cut by the domain edge differ, inside
/// the halo). Per level: aggregate (mean rock surface, mean debris, mean tracers), run `relax` on the coarse grid,
/// prolong the change back (losing cells scale their posts' debris and tracers by the cell ratio, exactly mass
/// conserving; gaining cells add the gain over their non-sink posts whose rock surface can hold debris, weighted
/// by the coarse new thickness interpolated at each post (smooth edges, exact per cell),
/// `accept` = rock slope <= repose + 5 deg, so coarse levels do not paint debris films onto faces; a gaining cell with
/// no such post spreads over its non-sink posts, and one with none loses the gain, counted). The last level must be
/// F = 1.
/// Influence radius = sum over levels of (2 N_L + 1) F_L dx.
#[allow(clippy::too_many_arguments)]
pub fn relax_multilevel(gx0: i64, gy0: i64, nx: usize, ny: usize, dx: f64, base: &[f64], accept: &[bool], tan_fine: Option<&[f64]>, t: &mut Vec<f64>, tracers: &mut [Vec<f64>], sink: &[u8], tan_r: f64, levels: &[(u32, u32)], beta: f64, threads: usize) -> RelaxOut {
    let mut lost = 0.0;
    let mut moved = 0.0;
    let mut its = 0;
    for &(f, n_it) in levels {
        if f <= 1 {
            let ro = relax_area(nx, ny, dx, base, None, tan_fine, t, tracers, sink, tan_r, n_it, beta, threads);
            lost += ro.lost_m3;
            moved += ro.moved_m3;
            its += n_it;
            continue;
        }
        let f = f as i64;
        let fd = |a: i64| a.div_euclid(f);
        let (c0x, c0y) = (fd(gx0), fd(gy0));
        let cnx = (fd(gx0 + nx as i64 - 1) - c0x + 1) as usize;
        let cny = (fd(gy0 + ny as i64 - 1) - c0y + 1) as usize;
        let cn = cnx * cny;
        let mut cnt = vec![0u32; cn];
        let mut cnt_free = vec![0u32; cn];
        let mut cnt_ns = vec![0u32; cn];
        let mut cb = vec![0.0f64; cn];
        let mut ct = vec![0.0f64; cn];
        let mut ctr: Vec<Vec<f64>> = tracers.iter().map(|_| vec![0.0f64; cn]).collect();
        let cell = |i: usize, j: usize| -> usize { ((fd(gy0 + j as i64) - c0y) as usize) * cnx + (fd(gx0 + i as i64) - c0x) as usize };
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let c = cell(i, j);
                cnt[c] += 1;
                if sink[k] == 0 && accept[k] { cnt_free[c] += 1; }
                if sink[k] == 0 { cnt_ns[c] += 1; }
                cb[c] += base[k];
                ct[c] += t[k];
                for (x, tr) in tracers.iter().enumerate() { ctr[x][c] += tr[k]; }
            }
        }
        for c in 0..cn {
            let w = cnt[c].max(1) as f64;
            cb[c] /= w;
            ct[c] /= w;
            for tr in ctr.iter_mut() { tr[c] /= w; }
        }
        let ct_old = ct.clone();
        let ctr_old = ctr.clone();
        let csink = vec![0u8; cn];
        let carea: Vec<f64> = cnt.iter().map(|&c| c as f64 / (f * f) as f64).collect();
        let rc = relax_area(cnx, cny, dx * f as f64, &cb, Some(&carea), None, &mut ct, &mut ctr, &csink, tan_r, n_it, beta, threads);
        moved += rc.moved_m3;
        its += n_it;
        // smooth, exactly conserving prolongation of gains: within a gaining cell the gain goes to its eligible posts
        // in proportion to the coarse new thickness interpolated (bilinear between coarse cell centres) at the post,
        // so a thin gain does not become a flat-topped 64 m sheet with a stepped edge
        let ff = f as f64;
        let ratio: Vec<f64> = (0..cn).map(|c| if ct_old[c] > 0.0 { ct[c] / ct_old[c] } else { 1.0 }).collect();
        // bilinear interpolation of a coarse field between coarse cell centres at fine post (i, j)
        let interp = |fld: &[f64], i: usize, j: usize| -> f64 {
            let u = ((gx0 + i as i64) as f64 - (ff - 1.0) / 2.0) / ff - c0x as f64;
            let v = ((gy0 + j as i64) as f64 - (ff - 1.0) / 2.0) / ff - c0y as f64;
            let u = u.clamp(0.0, (cnx - 1) as f64);
            let v = v.clamp(0.0, (cny - 1) as f64);
            let (iu, iv) = ((u.floor() as usize).min(cnx.saturating_sub(2)), (v.floor() as usize).min(cny.saturating_sub(2)));
            let (a, b) = (if cnx > 1 { u - iu as f64 } else { 0.0 }, if cny > 1 { v - iv as f64 } else { 0.0 });
            let at = |x: usize, y: usize| fld[y.min(cny - 1) * cnx + x.min(cnx - 1)];
            (1.0 - a) * (1.0 - b) * at(iu, iv) + a * (1.0 - b) * at(iu + 1, iv) + (1.0 - a) * b * at(iu, iv + 1) + a * b * at(iu + 1, iv + 1)
        };
        let wt = |i: usize, j: usize| -> f64 { interp(&ct, i, j).max(0.0) + 1.0e-9 };
        // losing cells: multiplicative with a smooth multiplier m = interpolated coarse ratio, renormalized per cell
        // so the cell keeps exactly its new coarse volume
        let mloss = |i: usize, j: usize| -> f64 { interp(&ratio, i, j).max(0.0) + 1.0e-9 };
        let eligible = |k: usize, c: usize| sink[k] == 0 && (accept[k] || cnt_free[c] == 0);
        let mut wsum = vec![0.0f64; cn];
        let mut lsum = vec![0.0f64; cn];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let c = cell(i, j);
                if ct[c] > ct_old[c] && eligible(k, c) { wsum[c] += wt(i, j); }
                if ct[c] <= ct_old[c] { lsum[c] += t[k] * mloss(i, j); }
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let c = cell(i, j);
                if ct[c] <= ct_old[c] {
                    let r = if lsum[c] > 0.0 { mloss(i, j) * ct[c] * cnt[c] as f64 / lsum[c] } else { 0.0 };
                    moved += (t[k] - t[k] * r) * dx * dx;
                    t[k] *= r;
                    for tr in tracers.iter_mut() { tr[k] *= r; }
                } else if eligible(k, c) {
                    let g = cnt[c] as f64 * wt(i, j) / wsum[c];
                    t[k] += (ct[c] - ct_old[c]) * g;
                    moved += (ct[c] - ct_old[c]) * g * dx * dx;
                    for (x, tr) in tracers.iter_mut().enumerate() { tr[k] = (tr[k] + (ctr[x][c] - ctr_old[x][c]) * g).max(0.0); }
                }
            }
        }
        for c in 0..cn {
            if ct[c] > ct_old[c] && cnt_ns[c] == 0 { lost += (ct[c] - ct_old[c]) * cnt[c] as f64 * dx * dx; }
        }
    }
    RelaxOut { lost_m3: lost, moved_m3: moved, iters: its }
}
