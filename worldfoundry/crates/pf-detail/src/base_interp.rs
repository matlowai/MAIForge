//! Opt-in C1, range-preserving reconstruction of published 32 m heights.
//!
//! Heights are geometry, not a second erosion solve. The four cell corners
//! bound every Bezier control ordinate, hence every reconstructed height.
//! Each vertex has ONE pair of tangents shared by all four incident patches.
//! Tangents start with the equal-spacing harmonic (Fritsch-Butland) limiter;
//! a single nonnegative scale then bounds all incident interior controls.
//! Mixed derivatives are zero. Shared values, first and mixed derivatives
//! make the bicubic patches C1 across both axes. No value clipping is used.
//! Edge curves are monotone between their endpoints. Interior values are
//! range-preserving; arbitrary oblique cuts need not be monotone.
//!
//! Two posts of real halo beyond queried cells are required; window-edge
//! clamping would make adjacent windows reconstruct different surfaces.
use crate::{input::T1bWindow, R};

pub const RULE: &str = "base_c1_v1: shared harmonic limited vertex tangents; zero mixed derivative; incident-cell Bezier convex-hull limiter; exact 32m posts; nearest-mm output; no new simulation time; no per-tile mean correction";

#[derive(Clone, Debug, serde::Serialize)]
pub struct Params { pub rule: &'static str }
impl Default for Params { fn default() -> Self { Self { rule: RULE } } }

pub struct C1Surface {
    gx0: i64,
    gy0: i64,
    n: usize,
    patches: Vec<[f64; 16]>,
}
fn tangent(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 || a.is_sign_positive() != b.is_sign_positive() { 0.0 }
    else { 2.0 * a * b / (a + b) }
}
fn bernstein(t: f64) -> [f64; 4] {
    let s = 1.0 - t;
    [s*s*s, 3.0*t*s*s, 3.0*t*t*s, t*t*t]
}
impl C1Surface {
    pub fn new(w: &T1bWindow) -> R<Self> {
        Self::from_grid(w.gx0, w.gy0, w.n, &w.height_m)
    }
    fn from_grid(gx0: i64, gy0: i64, n: usize, z: &[f64]) -> R<Self> {
        if n < 4 || z.len() != n*n || z.iter().any(|v| !v.is_finite()) {
            return Err("base_interp=c1 needs a finite square grid with real halo".into());
        }
        let mut g = vec![[0.0; 2]; n*n];
        for y in 1..n-1 { for x in 1..n-1 {
            let k = y*n+x; let h = z[k];
            let (dx,dy) = (tangent(h-z[k-1], z[k+1]-h), tangent(h-z[k-n], z[k+n]-h));
            let mut scale = 1.0_f64;
            for sy in [-1isize,1] { for sx in [-1isize,1] {
                let xn = (x as isize+sx) as usize; let yn = (y as isize+sy) as usize;
                let corners = [h,z[y*n+xn],z[yn*n+x],z[yn*n+xn]];
                let lo = corners.iter().copied().fold(f64::INFINITY,f64::min);
                let hi = corners.iter().copied().fold(f64::NEG_INFINITY,f64::max);
                // The edge controls are already bounded by the 1D limiter.
                // Bound the single interior control contributed by this vertex.
                let d = (sx as f64*dx+sy as f64*dy)/3.0;
                if d > 0.0 { scale = scale.min((hi-h)/d); }
                if d < 0.0 { scale = scale.min((lo-h)/d); }
            }}
            g[k] = [dx*scale,dy*scale];
        }}
        let m = n-3;
        let mut patches = vec![[0.0;16];m*m];
        for y in 1..n-2 { for x in 1..n-2 {
            let p = &mut patches[(y-1)*m+x-1];
            for cy in 0..2 { for cx in 0..2 {
                let k = (y+cy)*n+x+cx;
                let (bx,by) = (cx*3,cy*3);
                let (sx,sy) = (if cx==0 {1isize} else {-1},if cy==0 {1isize} else {-1});
                let (ix,iy) = ((bx as isize+sx) as usize,(by as isize+sy) as usize);
                p[by*4+bx] = z[k];
                p[by*4+ix] = z[k]+sx as f64*g[k][0]/3.0;
                p[iy*4+bx] = z[k]+sy as f64*g[k][1]/3.0;
                p[iy*4+ix] = z[k]+(sx as f64*g[k][0]+sy as f64*g[k][1])/3.0;
            }}
        }}
        Ok(Self {gx0:gx0+1,gy0:gy0+1,n:m,patches})
    }
    /// Coefficients in metres. Coordinates are GLOBAL 32 m cell addresses.
    pub fn patch(&self, gx: i64, gy: i64) -> &[f64;16] {
        let (x,y) = (gx-self.gx0,gy-self.gy0);
        assert!(x>=0 && y>=0 && x<self.n as i64 && y<self.n as i64,
            "base_interp=c1 outside real parent halo at ({gx},{gy}); refusing window-edge clamping");
        &self.patches[y as usize*self.n+x as usize]
    }
    pub fn eval_patch(p: &[f64;16], u: f64, v: f64) -> f64 {
        let a=bernstein(u); let b=bernstein(v);
        (0..4).map(|j|b[j]*(0..4).map(|i|a[i]*p[j*4+i]).sum::<f64>()).sum()
    }
    pub fn height_mm(&self, px: i64, py: i64) -> i32 {
        let p=self.patch(px.div_euclid(16),py.div_euclid(16));
        crate::to_mm(Self::eval_patch(p,px.rem_euclid(16) as f64/16.0,py.rem_euclid(16) as f64/16.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Numerical interpolation falsifiers, not geological-world fixtures.
    #[test]
    fn exact_bounded_and_c1_on_irregular_scalar_grid() {
        let n=19;
        let z:Vec<_>=(0..n*n).map(|k| (((k*7919+17)^(k*37))%20001) as f64/10.0-1000.0).collect();
        let s=C1Surface::from_grid(-9,-11,n,&z).unwrap();
        for y in 1..n-2 {for x in 1..n-2 {
            let p=s.patch(x as i64-9,y as i64-11);
            let c=[z[y*n+x],z[y*n+x+1],z[(y+1)*n+x],z[(y+1)*n+x+1]];
            let lo=c.iter().copied().fold(f64::INFINITY,f64::min);
            let hi=c.iter().copied().fold(f64::NEG_INFINITY,f64::max);
            for &v in p {assert!(v>=lo-1e-10 && v<=hi+1e-10);}
            for j in 0..=16 {for i in 0..=16 {
                let h=C1Surface::eval_patch(p,i as f64/16.0,j as f64/16.0);
                assert!(h>=lo-1e-10 && h<=hi+1e-10);
            }}
            assert_eq!(s.height_mm((x as i64-9)*16,(y as i64-11)*16),crate::to_mm(z[y*n+x]));
            if x<n-3 {let q=s.patch(x as i64-8,y as i64-11);for j in 0..4 {
                assert_eq!(p[j*4+3],q[j*4]);
                assert!(((p[j*4+3]-p[j*4+2])-(q[j*4+1]-q[j*4])).abs()<1e-10);
            }}
            if y<n-3 {let q=s.patch(x as i64-9,y as i64-10);for i in 0..4 {
                assert_eq!(p[12+i],q[i]);
                assert!(((p[12+i]-p[8+i])-(q[4+i]-q[i])).abs()<1e-10);
            }}
        }}
    }
    #[test]
    fn planar_reproduction_and_monotone_edges() {
        let n=10;let z:Vec<_>=(0..n*n).map(|k|3.0*(k%n) as f64-7.0*(k/n) as f64).collect();
        let s=C1Surface::from_grid(0,0,n,&z).unwrap();
        for j in 16..112 {for i in 16..112 {
            assert!((s.height_mm(i,j) as f64/1000.0-(3.0*i as f64-7.0*j as f64)/16.0).abs()<=0.000501);
        }}
        let z:Vec<_>=(0..n*n).map(|k|((k%n).pow(3)+(k/n).pow(2)) as f64).collect();
        let s=C1Surface::from_grid(0,0,n,&z).unwrap();
        for gy in 1..8 {for px in 16..111 {assert!(s.height_mm(px+1,gy*16)>=s.height_mm(px,gy*16));}}
    }
    #[test]
    fn overlapping_real_support_is_independent_of_window_origin() {
        let n=15;let z:Vec<_>=(0..n*n).map(|k|(k as f64*0.31).sin()*20.0).collect();
        let a=C1Surface::from_grid(-7,-7,n,&z).unwrap();
        let z2:Vec<_>=(0..10).flat_map(|y|(0..10).map({let z=&z;move |x|z[(y+3)*n+x+2]})).collect();
        let b=C1Surface::from_grid(-5,-4,10,&z2).unwrap();
        for py in -48..64 {for px in -64..48 {assert_eq!(a.height_mm(px,py),b.height_mm(px,py));}}
    }
    #[test]
    #[should_panic(expected="refusing window-edge clamping")]
    fn missing_halo_is_not_silently_clamped() {
        let s=C1Surface::from_grid(0,0,4,&[0.0;16]).unwrap();s.height_mm(0,0);
    }
    #[test]
    fn key_off_omits_recipe_field_and_bad_values_fail() {
        for mut p in [crate::DetailParamsV1::v1(),crate::DetailParamsV1::v2(),crate::DetailParamsV1::v3()] {
            let off=serde_json::to_value(&p).unwrap();assert!(off.get("base_interp").is_none());
            p.set("base_interp","c1").unwrap();assert!(serde_json::to_value(&p).unwrap().get("base_interp").is_some());
            p.set("base_interp","off").unwrap();assert_eq!(serde_json::to_value(&p).unwrap(),off);
            assert!(p.set("base_interp","catmull").is_err());
        }
    }
}
