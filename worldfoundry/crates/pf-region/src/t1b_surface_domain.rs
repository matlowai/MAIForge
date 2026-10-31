//! Canonical detail-domain inputs for the LS routing owner.
//!
//! A crop edge is never a terminal. A component is closed only when every legal
//! neighbour of every active post is present, terminating at declared native
//! boundaries. Inputs are sorted by global (gy,gx), independent of arrival order.
//! This conservative closure is sufficient, not necessary: an open component
//! yields NeedsDomain rather than an approximation or accepted output.
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub const VERSION: &str = "ls_detail_component_v1";
pub const D8: [(i64, i64); 8] = [(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)];

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PostId { pub gy: i64, pub gx: i64 }
impl PostId {
    pub fn neighbour(self, k: usize) -> Self { Self { gy: self.gy+D8[k].1, gx: self.gx+D8[k].0 } }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct StagePin {
    pub generation: String,
    pub chart: String,
    pub recipe: String,
    pub stage: String,
    pub predecessor: String,
    pub ground: String,
    pub boundary: String,
    pub sources: String,
    pub spacing_mm: i64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Policy {
    PinnedConstructionV1,
    /// Historical explicitly basin-confined detail contract.
    NativeDetailV1,
    /// Actual physical fixed boundaries and declared legal/conduit edges govern
    /// the mutable solve. A coarse catchment label is provenance, not a wall.
    /// Physical MFD retains every permitted downhill D8 flux, including crossing
    /// diagonal fractions; the legacy noncrossing line constraint is not applied.
    NativePhysicalV1,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mode { SingleFlowV1, MultipleFlowSlope11V1 }
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Boundary {
    pub id: String,
    pub basin: u32,
    /// Only CHANNEL, GUIDE, MARINE and LAKE. RING is never supported.
    pub kind: u8,
    /// Native network attachment for CHANNEL/GUIDE; absent for native water.
    pub attachment: Option<u32>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Post {
    pub id: PostId,
    pub ground_m: f64,
    pub support_m2: f64,
    /// Own runoff-effective source, never accumulated area.
    pub source_m2: f64,
    pub external_m2: f64,
    pub basin: u32,
    /// Outgoing legal D8 edges supplied by the owner, not a visual guide mask.
    pub legal: u8,
    pub locked_receiver: Option<PostId>,
    pub boundary: Option<Boundary>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tile {
    pub pin: StagePin,
    pub policy: Policy,
    pub posts: Vec<Post>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum DomainError {
    NeedsDomain { stage: String, frontier: Vec<PostId> },
    Stale { field: String },
    DuplicateSource { post: PostId },
    Invalid { reason: String },
    Resource { posts: usize, limit: usize },
}
impl std::fmt::Display for DomainError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "{}", serde_json::to_string(self).unwrap()) }
}
impl std::error::Error for DomainError {}
pub type Result<T> = std::result::Result<T, DomainError>;
pub fn invalid(reason: impl Into<String>) -> DomainError { DomainError::Invalid { reason: reason.into() } }

/// Private fields ensure routing can only consume validated, closed inputs.
#[derive(Clone, Debug)]
pub struct Domain {
    pub(crate) pin: StagePin,
    pub(crate) policy: Policy,
    pub(crate) posts: Vec<Post>,
    pub(crate) index: BTreeMap<PostId, usize>,
    pub(crate) digest: String,
}
impl Domain {
    pub fn pin(&self) -> &StagePin { &self.pin }
    pub fn posts(&self) -> &[Post] { &self.posts }
    pub fn digest(&self) -> &str { &self.digest }
    pub fn index(&self, id: PostId) -> Option<usize> { self.index.get(&id).copied() }
    pub fn assemble(pin: &StagePin, policy: Policy, tiles: Vec<Tile>, max_posts: usize) -> Result<Self> {
        if pin.spacing_mm <= 0 || [&pin.generation,&pin.chart,&pin.recipe,&pin.stage,&pin.predecessor,&pin.ground,&pin.boundary,&pin.sources].iter().any(|s| s.is_empty()) {
            return Err(invalid("complete content pins and positive spacing required"));
        }
        let mut by_id = BTreeMap::new();
        for tile in tiles {
            if tile.pin != *pin || tile.policy != policy { return Err(DomainError::Stale { field: "tile stage/recipe/predecessor/boundary/source/policy".into() }); }
            for p in tile.posts {
                let id=p.id;
                if by_id.insert(id,p).is_some() { return Err(DomainError::DuplicateSource { post:id }); }
                if by_id.len()>max_posts { return Err(DomainError::Resource { posts:by_id.len(),limit:max_posts }); }
            }
        }
        if by_id.is_empty() { return Err(invalid("empty domain")); }
        let posts:Vec<Post>=by_id.into_values().collect();
        let index:BTreeMap<PostId,usize>=posts.iter().enumerate().map(|(i,p)|(p.id,i)).collect();
        let mut missing=BTreeSet::new();
        let mut boundary_defs=BTreeMap::new();
        for p in &posts {
            if !p.ground_m.is_finite() || !p.support_m2.is_finite() || p.support_m2<=0.0 || !p.source_m2.is_finite() || p.source_m2<0.0 || !p.external_m2.is_finite() || p.external_m2<0.0 {
                return Err(invalid(format!("nonfinite height or invalid source/support at {:?}",p.id)));
            }
            if let Some(b)=&p.boundary {
                use crate::t1b::{CHANNEL,GUIDE,MARINE,LAKE};
                if b.id.is_empty() || ![CHANNEL,GUIDE,MARINE,LAKE].contains(&b.kind) || b.basin!=p.basin || p.locked_receiver.is_some() || p.source_m2!=0.0 || p.external_m2!=0.0 {
                    return Err(invalid(format!("unsupported terminal or source at boundary {:?}",p.id)));
                }
                if [CHANNEL,GUIDE].contains(&b.kind) && b.attachment.is_none() { return Err(invalid("channel boundary lacks native attachment")); }
                if let Some(old)=boundary_defs.insert(b.id.clone(),b.clone()) { if old!=*b { return Err(invalid("conflicting boundary identity")); } }
                continue;
            }
            if p.legal==0 { return Err(invalid(format!("active post has no legal edge {:?}",p.id))); }
            if policy==Policy::PinnedConstructionV1 && (p.legal!=255 || p.locked_receiver.is_some()) { return Err(invalid("construction policy requires original unrestricted D8 semantics")); }
            if let Some(q)=p.locked_receiver {
                if !(0..8).any(|k|p.legal&(1<<k)!=0 && p.id.neighbour(k)==q) { return Err(invalid("locked receiver is not a legal neighbour")); }
            }
            for k in 0..8 { if p.legal&(1<<k)!=0 {
                let q=p.id.neighbour(k);
                if !index.contains_key(&q) { missing.insert(q); }
            } }
        }
        if !missing.is_empty() { return Err(DomainError::NeedsDomain {stage:pin.stage.clone(),frontier:missing.into_iter().collect()}); }
        let (ground,sources,boundary)=component_digests(&posts)?;
        for (field,actual,expected) in [("ground",ground,&pin.ground),("sources",sources,&pin.sources),("boundary",boundary,&pin.boundary)] {
            if actual!=*expected {return Err(DomainError::Stale {field:field.into()});}
        }
        let digest=content_digest(&(VERSION,pin,policy,&posts))?;
        Ok(Self { pin:pin.clone(),policy,posts,index,digest })
    }
    pub(crate) fn edge(&self, i:usize, k:usize) -> Option<usize> {
        let p=&self.posts[i];
        if p.boundary.is_some() || p.legal&(1<<k)==0 { return None; }
        let id=p.id.neighbour(k);
        if p.locked_receiver.is_some_and(|q|q!=id) { return None; }
        let j=*self.index.get(&id)?;
        if self.policy==Policy::NativeDetailV1 && p.basin!=self.posts[j].basin { return None; }
        Some(j)
    }
}
pub fn content_digest<T: Serialize>(value:&T)->Result<String> {
    struct HashWriter(blake3::Hasher);
    impl std::io::Write for HashWriter {
        fn write(&mut self,b:&[u8])->std::io::Result<usize>{self.0.update(b);Ok(b.len())}
        fn flush(&mut self)->std::io::Result<()>{Ok(())}
    }
    let mut out=HashWriter(blake3::Hasher::new());
    serde_json::to_writer(&mut out,value).map_err(|e|invalid(e.to_string()))?;
    Ok(out.0.finalize().to_hex().to_string())
}
/// Path is a locator only; the expected digest is mandatory and checked first.
pub fn read_tile(path:&Path, expected:&str)->Result<Tile> {
    let bytes=std::fs::read(path).map_err(|e|invalid(e.to_string()))?;
    if blake3::hash(&bytes).to_hex().as_str()!=expected { return Err(DomainError::Stale {field:"tile bytes".into()}); }
    serde_json::from_slice(&bytes).map_err(|e|invalid(e.to_string()))
}

/// Exact, read-only stage capture. This is diagnostic source evidence, not a
/// canonical certificate: prior fill taint and evaluation-ring dependency remain
/// explicit. Never reconstruct pre-fill ground from a filled published window.
pub fn capture_stage(lat:&crate::t1b_rivernet::EvalLat,z:&[f64],outlet:&[bool],sm:&crate::t1b_seam::SeamTaint,stage:&str)->std::result::Result<(),String> {
    let Some(root)=std::env::var_os("VZ_SURFACE_CAPTURE") else {return Ok(())};
    let root=std::path::PathBuf::from(root);
    let name=format!("{stage}_{}_{}_{}",lat.gx0,lat.gy0,lat.n);
    let dir=root.join(name);
    std::fs::create_dir_all(&root).map_err(|e|e.to_string())?;
    std::fs::create_dir(&dir).map_err(|e|format!("surface capture refuses overwrite: {e}"))?;
    let mut fields=BTreeMap::new();
    for (name,bytes) in [("ground.f64",z.iter().flat_map(|v|v.to_le_bytes()).collect::<Vec<_>>()),("outlet.u8",outlet.iter().map(|v|*v as u8).collect()),("prior_taint.u8",sm.taint.iter().map(|v|*v as u8).collect())] {
        let digest=blake3::hash(&bytes).to_hex().to_string();
        std::fs::write(dir.join(name),bytes).map_err(|e|e.to_string())?;
        fields.insert(name,digest);
    }
    let meta=serde_json::json!({"schema":"ls-stage-capture/v1","stage":stage,"gx0":lat.gx0,"gy0":lat.gy0,"n":lat.n,"spacing_mm":lat.dx_mm,"margin_posts":sm.margin_posts,"fields":fields,"accepted":false,"note":"Uncertified exact pre-fill evidence. Outlet mask alone does not establish native terminal identity or source closure."});
    std::fs::write(dir.join("capture.json"),serde_json::to_vec_pretty(&meta).unwrap()).map_err(|e|e.to_string())
}

/// Exact boundary classification beside the pre-fill snapshot. The mask is not
/// changed. Native channel attachments still need the upstream network proof.
pub fn capture_fixed(lat:&crate::t1b_rivernet::EvalLat,fixed:&[u8],stage:&str)->std::result::Result<(),String> {
    let Some(root)=std::env::var_os("VZ_SURFACE_CAPTURE") else {return Ok(())};
    let dir=std::path::PathBuf::from(root).join(format!("{stage}_{}_{}_{}",lat.gx0,lat.gy0,lat.n));
    let path=dir.join("fixed.u8");
    use std::io::Write;
    std::fs::OpenOptions::new().write(true).create_new(true).open(path).and_then(|mut f|f.write_all(fixed)).map_err(|e|e.to_string())?;
    let file=dir.join("capture.json");
    let mut meta:serde_json::Value=serde_json::from_slice(&std::fs::read(&file).map_err(|e|e.to_string())?).map_err(|e|e.to_string())?;
    meta["fields"]["fixed.u8"]=blake3::hash(fixed).to_hex().to_string().into();
    std::fs::write(file,serde_json::to_vec_pretty(&meta).unwrap()).map_err(|e|e.to_string())
}

/// Read-only source and attachment evidence for stage qualification. Every field
/// retains global addresses and content digests; these are not accepted tiles.
pub fn capture_native(ctx:&crate::t1b_rivernet::RnContext<'_>,lp:&crate::t1b_rivernet::LatPrep<'_>,s:&crate::t1b_rivernet::Surface)->std::result::Result<(),String> {
    let Some(root)=std::env::var_os("VZ_SURFACE_CAPTURE") else {return Ok(())};
    let lat=&lp.lat;
    let dir=std::path::PathBuf::from(root).join(format!("native_{}_{}_{}",lat.gx0,lat.gy0,lat.n));
    std::fs::create_dir_all(dir.parent().unwrap()).map_err(|e|e.to_string())?;
    std::fs::create_dir(&dir).map_err(|e|e.to_string())?;
    let mut fields=BTreeMap::new();
    let mut write=|name:&str,b:Vec<u8>|->std::result::Result<(),String>{
        fields.insert(name.to_string(),blake3::hash(&b).to_hex().to_string());
        std::fs::write(dir.join(name),b).map_err(|e|e.to_string())
    };
    write("fixed.u8",s.fixed.clone())?;
    write("chan_link.u32",s.chan_link.iter().flat_map(|v|v.to_le_bytes()).collect())?;
    write("cell.u32",lp.cell.iter().flat_map(|v|v.to_le_bytes()).collect())?;
    let mut basins=BTreeMap::new();
    for &c in &lp.cell {if !basins.contains_key(&c) {basins.insert(c,ctx.t0.final_basin(c)?);}}
    write("basin.u32",lp.cell.iter().flat_map(|c|basins[c].to_le_bytes()).collect())?;
    write("source.f64",(0..lat.nodes()).flat_map(|i|(lat.dx()*lat.dx()*if lp.rw.is_empty(){1.0}else{lp.rw[i] as f64}).to_le_bytes()).collect())?;
    let meta=serde_json::json!({"schema":"ls-native-capture/v1","gx0":lat.gx0,"gy0":lat.gy0,"n":lat.n,"spacing_mm":lat.dx_mm,"generation":ctx.rn_recipe_digest,"layer":ctx.layer_digest,"fields":fields,"accepted":false,"note":"Local source and attachment evidence; boundary/confinement qualification and canonical upstream stage closure still required."});
    std::fs::write(dir.join("capture.json"),serde_json::to_vec_pretty(&meta).unwrap()).map_err(|e|e.to_string())
}

#[derive(Deserialize)]
struct TileInput { path:String, digest:String }
#[derive(Deserialize)]
struct BuildRequest { pin:StagePin, policy:Policy, mode:Mode, max_posts:usize, tiles:Vec<TileInput> }
/// Build a candidate owner artifact from content-bound, single-owner stage tiles.
/// This does not accept it for game publication: F/seam gates remain external.
pub fn build_from_request(request:&Path,out:&Path)->Result<()> {
    let bytes=std::fs::read(request).map_err(|e|invalid(e.to_string()))?;
    let req:BuildRequest=serde_json::from_slice(&bytes).map_err(|e|invalid(e.to_string()))?;
    if req.max_posts==0 || req.max_posts>5_000_000 {return Err(invalid("surface-build max_posts must be 1..5000000 under the 8 GB lane cap"));}
    std::fs::create_dir(out).map_err(|e|invalid(format!("surface-build refuses overwrite: {e}")))?;
    let result=(||{
        let mut tiles=Vec::new();
        let mut count=0usize;
        for entry in req.tiles {
            let path=request.parent().unwrap_or(Path::new(".")).join(entry.path);
            let tile=read_tile(&path,&entry.digest)?;
            count=count.checked_add(tile.posts.len()).ok_or_else(||invalid("post count overflow"))?;
            if count>req.max_posts {return Err(DomainError::Resource {posts:count,limit:req.max_posts});}
            tiles.push(tile);
        }
        let domain=Domain::assemble(&req.pin,req.policy,tiles,req.max_posts)?;
        let route=crate::t1b_ls::route_detail_surface(&domain,req.mode)?;
        std::fs::write(out.join("domain.json"),serde_json::to_vec(&Tile {pin:domain.pin.clone(),policy:domain.policy,posts:domain.posts.clone()}).map_err(|e|invalid(e.to_string()))?).map_err(|e|invalid(e.to_string()))?;
        std::fs::write(out.join("route.json"),serde_json::to_vec(&route).map_err(|e|invalid(e.to_string()))?).map_err(|e|invalid(e.to_string()))?;
        let cert=serde_json::json!({"schema":VERSION,"pin":req.pin,"policy":req.policy,"domain_digest":domain.digest,"route_digest":route.digest,"posts":count,"closed":true,"accepted_for_publication":false,"unmet_gates":["upstream-stage qualification","native window integration","literal F1-F7","full seam matrix","morphology","Unreal"]});
        std::fs::write(out.join("certificate.json"),serde_json::to_vec_pretty(&cert).unwrap()).map_err(|e|invalid(e.to_string()))?;
        Ok(())
    })();
    if let Err(e)=&result {
        std::fs::write(out.join("failure.json"),serde_json::to_vec_pretty(e).unwrap()).map_err(|io|invalid(io.to_string()))?;
    }
    result
}

/// Digests of exact component-owned fields, independent of partition. Canonical
/// source producers pin these once; consumers compare instead of repinning.
pub fn component_digests(posts:&[Post])->Result<(String,String,String)> {
    let mut ordered:Vec<&Post>=posts.iter().collect();ordered.sort_by_key(|p|p.id);
    let ground:Vec<_>=ordered.iter().map(|p|(p.id,p.ground_m)).collect();
    let sources:Vec<_>=ordered.iter().map(|p|(p.id,p.support_m2,p.source_m2,p.external_m2)).collect();
    let boundary:Vec<_>=ordered.iter().map(|p|(p.id,p.basin,p.legal,p.locked_receiver,&p.boundary)).collect();
    Ok((content_digest(&ground)?,content_digest(&sources)?,content_digest(&boundary)?))
}
