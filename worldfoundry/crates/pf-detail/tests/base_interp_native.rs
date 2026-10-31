//! Real-data numerical certificate. Set VZ_C1_INPUT to a small manifest of
//! published windows; no implicit dependency on a developer's local paths.
use pf_detail::{base_interp::C1Surface,input::T1bWindow};
use serde_json::{Value,json};
use std::path::Path;
fn latlon(h:&Value,x:f64,y:f64)->(f64,f64) {
    let e=&h["prep"]["net"]["evolved_surface"];let c=&e["chart"];
    let r=c["datum_radius_mm"].as_f64().unwrap()/1000.;
    let a=x/(2.*r);let b=y/(2.*r);let q=a*a+b*b;
    let (ve,vn,vu)=(2.*a/(1.+q),2.*b/(1.+q),(1.-q)/(1.+q));
    let la0=(c["origin_lat_e7"].as_f64().unwrap()/1e7).to_radians();
    let lo0=(c["origin_lon_e7"].as_f64().unwrap()/1e7).to_radians();
    let px=la0.cos()*vu-la0.sin()*vn;let pz=la0.sin()*vu+la0.cos()*vn;
    let la=pz.atan2(px.hypot(ve)).to_degrees();let lo=(lo0+ve.atan2(px)).to_degrees();
    (la,lo)
}
#[test]
fn native_bounds_posts_derivative_seams_whole_cell_and_window_overlap() {
    let Ok(path)=std::env::var("VZ_C1_INPUT")else{eprintln!("SKIP native C1 certificate: set VZ_C1_INPUT");return};
    let m:Value=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    let dirs=m["windows"].as_array().unwrap();
    let mut reports=Vec::new();
    for (wi,d) in dirs.iter().enumerate() {
        let dir=Path::new(d.as_str().unwrap());let w=T1bWindow::load(dir).unwrap();let s=C1Surface::new(&w).unwrap();
        let mut anchors=0u64;let mut controls=0u64;let mut seams=0u64;let mut max_error=0.0_f64;
        for y in 1..w.n-2 {for x in 1..w.n-2 {
            let (gx,gy)=(w.gx0+x as i64,w.gy0+y as i64);let p=s.patch(gx,gy);
            let c=[w.h(gx,gy),w.h(gx+1,gy),w.h(gx,gy+1),w.h(gx+1,gy+1)];
            let lo=c.iter().copied().fold(f64::INFINITY,f64::min);let hi=c.iter().copied().fold(f64::NEG_INFINITY,f64::max);
            for &v in p {assert!(v>=lo-1e-9 && v<=hi+1e-9,"overshoot {gx} {gy} {v} [{lo},{hi}]");controls+=1;}
            assert_eq!(s.height_mm(gx*16,gy*16),pf_detail::to_mm(w.h(gx,gy)));anchors+=1;
            if x<w.n-3 {let q=s.patch(gx+1,gy);for j in 0..4 {assert_eq!(p[j*4+3],q[j*4]);let err=((p[j*4+3]-p[j*4+2])-(q[j*4+1]-q[j*4])).abs()*3./32.;max_error=max_error.max(err);assert!(err<1e-10);seams+=1;}}
            if y<w.n-3 {let q=s.patch(gx,gy+1);for i in 0..4 {assert_eq!(p[12+i],q[i]);let err=((p[12+i]-p[8+i])-(q[4+i]-q[i])).abs()*3./32.;max_error=max_error.max(err);assert!(err<1e-10);seams+=1;}}
        }}
        let mut result=json!({"window":d,"digest":w.digest,"exact_parent_posts":anchors,"bounded_controls":controls,"derivative_seam_controls":seams,"derivative_error_max_m_per_m":max_error});
        if Some(wi as u64)==m["whole_window_index"].as_u64() {
            let doc:Value=serde_json::from_slice(&std::fs::read(dir.join("window.json")).unwrap()).unwrap();let h=&doc["header"];
            let target=m["whole_cell"].as_u64().unwrap() as u32;
            let mut xa=i64::MAX;let mut ya=i64::MAX;let mut xb=i64::MIN;let mut yb=i64::MIN;
            for y in 0..w.n {for x in 0..w.n {if w.t0_cell[y*w.n+x]==target {xa=xa.min(x as i64);xb=xb.max(x as i64);ya=ya.min(y as i64);yb=yb.max(y as i64);}}}
            assert!(xa>=3 && ya>=3 && xb<w.n as i64-4 && yb<w.n as i64-4,"whole target cell not enclosed in real halo");
            let ev=&h["prep"]["net"]["evolved_surface"];
            let col=ev["columns"].as_array().unwrap().iter().find(|c|c["cell"].as_u64()==Some(target as u64)).unwrap();
            let reference=col["mean_m"].as_f64().unwrap();
            let nr=col["shape"][0].as_u64().unwrap() as u32;let nc=col["shape"][1].as_u64().unwrap() as u32;
            let r=ev["chart"]["datum_radius_mm"].as_f64().unwrap()/1000.;
            let lat0=-90.+(target/nc) as f64*180./nr as f64;let lat1=lat0+180./nr as f64;
            let lonc=-180.+((target%nc) as f64+0.5)*360./nc as f64;let lonhalf=180./nc as f64;
            let area=r*r*(360./nc as f64).to_radians()*(lat1.to_radians().sin()-lat0.to_radians().sin());
            let mut lo=[0.0;2];let mut hi=[0.0;2];let mut area_lo=0.;let mut area_hi=0.;let mut inside=0u64;let mut boundary=0u64;
            for gy in w.gy0+ya-2..w.gy0+yb+2 {for gx in w.gx0+xa-2..w.gx0+xb+2 {
                let (x,y)=((gx as f64+0.5)*32.,(gy as f64+0.5)*32.);
                let (la,lon)=latlon(h,x,y);let dl=(lon-lonc+180.).rem_euclid(360.)-180.;
                // Stereographic metric = 1/(1+rho^2/(4R^2)) <= 1.
                // The entire square is inside this geodesic ball. Latitude
                // is 1/R Lipschitz; longitude is 1/(R cos(latitude)). These
                // conservative bounds classify whole patches without treating
                // a centre sample or four corners as exact cell membership.
                let dlat=(16.*2.0_f64.sqrt()/r).to_degrees()+1e-10;
                let dlon=dlat/(la.abs()+dlat).to_radians().cos();
                if la+dlat<lat0 || la-dlat>lat1 || dl+dlon < -lonhalf || dl-dlon>lonhalf {continue;}
                let complete=la-dlat>=lat0 && la+dlat<=lat1 && dl-dlon>=-lonhalf && dl+dlon<=lonhalf;
                let p=s.patch(gx,gy);let c=[p[0],p[3],p[12],p[15]];
                let mn=c.iter().copied().fold(f64::INFINITY,f64::min)-reference;
                let mx=c.iter().copied().fold(f64::NEG_INFINITY,f64::max)-reference;
                let rho_min=(x.abs()-16.).max(0.).powi(2)+(y.abs()-16.).max(0.).powi(2);
                let rho_max=(x.abs()+16.).powi(2)+(y.abs()+16.).powi(2);
                let wmin=1./(1.+rho_max/(4.*r*r)).powi(2);let wmax=1./(1.+rho_min/(4.*r*r)).powi(2);
                if complete {
                    inside+=1;area_lo+=1024.*wmin;area_hi+=1024.*wmax;
                    // Bernstein basis integrals are exactly 1/4 in each axis.
                    // Integrate the bicubic analytically, bound all area-Jacobian
                    // variation, and reserve 0.5 mm for integer publication.
                    let mean=[c.iter().sum::<f64>()/4.-reference,p.iter().sum::<f64>()/16.-reference];
                    let error=mn.abs().max(mx.abs())*(wmax-wmin)*512.;
                    for k in 0..2 {let mass=mean[k]*(wmin+wmax)*512.;lo[k]+=mass-error;hi[k]+=mass+error;}
                } else {
                    boundary+=1;area_hi+=1024.*wmax;
                    // Unknown exact area fraction is enclosed, never guessed.
                    for k in 0..2 {lo[k]+=mn.min(0.)*1024.*wmax;hi[k]+=mx.max(0.)*1024.*wmax;}
                }
            }}
            assert!(inside>0 && area>=area_lo && area<=area_hi,"whole-cell coverage certificate failed");
            let intervals:[[f64;2];2]=std::array::from_fn(|k|[reference+lo[k]/area-0.0005,reference+hi[k]/area+0.0005]);
            for interval in intervals {assert!(interval[0]>=reference-50. && interval[1]<=reference+50.,"F6 enclosure outside tolerance: {interval:?}");}
            result["whole_cell_f6"]=json!({"cell":target,"whole_patches":inside,"boundary_patches_enclosed":boundary,"native_area_m2":area,"covered_area_interval_m2":[area_lo,area_hi],"mean_before_interval_m":intervals[0],"mean_after_interval_m":intervals[1],"t0_mean_m":reference,"tolerance_m":50.,"pass":true,"method":"analytic bicubic integrals; rigorous Jacobian and geographic-boundary enclosures; plus 0.5mm publication quantization; no partial-window mean substituted"});
        }
        reports.push(result);
    }
    // Separate adjacent real windows must use identical shared support, not
    // merely two calls through the same context.
    let a=T1bWindow::load(Path::new(dirs[1].as_str().unwrap())).unwrap();let b=T1bWindow::load(Path::new(dirs[2].as_str().unwrap())).unwrap();
    let sa=C1Surface::new(&a).unwrap();let sb=C1Surface::new(&b).unwrap();
    let x0=(a.gx0+1).max(b.gx0+1)*16;let x1=(a.gx0+a.n as i64-2).min(b.gx0+b.n as i64-2)*16;
    let y0=(a.gy0+1).max(b.gy0+1)*16;let y1=(a.gy0+a.n as i64-2).min(b.gy0+b.n as i64-2)*16;
    let mut shared=0u64;
    for py in y0..y1 {for px in x0..x1 {assert_eq!(sa.height_mm(px,py),sb.height_mm(px,py),"native cross-window mismatch {px} {py}");shared+=1;}}
    assert!(shared>0);
    let report=json!({"windows":reports,"cross_window_2m_posts_exact":shared,"quantization":"nearest mm; C1 certificate is for continuous reconstruction before quantization"});
    println!("{}",serde_json::to_string_pretty(&report).unwrap());
    if let Some(out)=m["report"].as_str(){std::fs::write(out,serde_json::to_vec_pretty(&report).unwrap()).unwrap();}
}
