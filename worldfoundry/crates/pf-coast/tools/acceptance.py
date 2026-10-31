#!/usr/bin/env python3
"""Actual-input coast overlap/plots and GLO30 morphology proxy. No mocked gates.
Requires numpy, Pillow and rasterio (only for the real DEM comparison).
"""
import argparse, hashlib, json, math, pathlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
P=pathlib.Path
NAMES=['none','rocky cliff','rocky platform','gravel beach','sand beach','tidal flat','salt marsh potential','delta front potential']
COLORS=['#ffffff','#d95750','#925ca7','#a87d4d','#f4cf58','#57c8d2','#57ad65','#ec8ce0']
DT={'u8':'u1','u16':'<u2','f32':'<f4'}
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def font(n):return ImageFont.truetype(FONT,n)
def read(d,f,n):return np.fromfile(d/f['file'],DT[f['dtype']]).reshape(n,n)
def slope(z,dx,dy):
 s=np.zeros_like(z,dtype=np.float32)
 for yy in [-1,0,1]:
  for xx in [-1,0,1]:
   if yy==xx==0:continue
   q=np.roll(np.roll(z,yy,0),xx,1);s=np.maximum(s,np.abs(z-q)/math.hypot(dx*xx,dy*yy))
 s[[0,-1],:]=np.nan;s[:,[0,-1]]=np.nan
 return s

def paint(z,classes,path,title,sub,spacing=32,dy=None,real=False,fixed=None):
 dy=spacing if dy is None else dy
 gy,gx=np.gradient(z.astype(float),dy,spacing);shade=np.clip((0.6-0.5*gx+0.5*gy)/np.sqrt(1+gx*gx+gy*gy),0,1)
 gray=(65+shade*180).astype(np.uint8);rgb=np.stack([gray,gray,gray],axis=-1)
 if fixed is not None:
  rgb[fixed==2]=[85,135,175];rgb[(fixed==1)|(fixed==5)]=[90,165,205];rgb[fixed==4]=[90,145,190]
  low=(fixed==0)&(z>=-1)&(z<=1);rgb[low]=(0.5*rgb[low]+0.5*np.array([215,190,140])).astype(np.uint8)
 for i,c in enumerate(COLORS[1:],1):
  mask=classes==i;col=np.array(tuple(bytes.fromhex(c[1:])),dtype=float);rgb[mask]=(0.8*col+0.2*rgb[mask]).astype(np.uint8)
 aspect=z.shape[1]*spacing/(z.shape[0]*dy);size=(1100,round(1100/aspect)) if aspect>=1 else (round(1100*aspect),1100)
 mapim=Image.fromarray(np.flipud(rgb)).resize(size,Image.Resampling.NEAREST)
 im=Image.new('RGB',(1440,1250),'white');im.paste(mapim,(30,110));d=ImageDraw.Draw(im)
 d.text((30,18),title,font=font(27),fill='black');d.text((30,62),sub,font=font(15),fill='black')
 for i in ([1,2,4] if real else range(1,8)):
  y=130+i*55;d.rectangle((1150,y,1172,y+22),fill=COLORS[i]);d.text((1180,y),({1:'slope >=30 deg',2:'slope 10-30 deg',4:'slope <10 deg'}[i] if real else NAMES[i]),font=font(15),fill='black')
 d.text((1150,590),('Blue = native water' if fixed is not None else ''),font=font(15),fill='black');d.text((1150,620),('Tan = low native ground' if fixed is not None else ''),font=font(15),fill='black')
 d.text((1150,650),'North up',font=font(18),fill='black');d.text((1150,690),f'{z.shape[1]*spacing/1000:.1f} km across',font=font(16),fill='black')
 d.text((30,1220),('GLO30 surface model. Zero-elevation sea proxy. No shore-type ground truth or sediment data.' if real else '32 m posts. Narrow bands remain narrow. Potential classes do not assert sediment deposits or vegetation.'),font=font(15),fill='black')
 im.save(path)

def coast_site(inp,out,tag,notes):
 wins={};seams=[];manifest=[];z=np.zeros((1536,1536),np.float32);cls=np.zeros_like(z,np.uint8);sl=np.zeros_like(z);shore=np.zeros_like(z);fixed=np.zeros_like(z,np.uint8);flagcounts={}
 for d in sorted(out.glob('w_*')):
  h=json.load(open(d/'coast.json'))['header'];l=h['lattice'];wins[l['rx'],l['ry']]=(d,h)
  native=json.load(open(inp/d.name/'window.json'))['header'];nl=native['lattice'];assert native['prep']['net']['margin_mm']==8192000
  fields={f['name']:read(d,f,l['n']) for f in h['fields']}
  for f in h['fields']:
   b=(d/f['file']).read_bytes();manifest.append({'file':str(d/f['file']),'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
  for k in ['coast_width_m','sediment_supply_index','wave_height_m','shore_slope']:assert np.isfinite(fields[k]).all()
  assert (fields['sediment_supply_index']>=0).all() and (fields['sediment_supply_index']<=1).all()
  assert (fields['shore_class']<8).all()
  yy=(l['ry']+1)*512;xx=(l['rx']+1)*512;core=np.s_[16:528,16:528]
  cls[yy:yy+512,xx:xx+512]=fields['shore_class'][core];sl[yy:yy+512,xx:xx+512]=fields['shore_slope'][core]
  nh=nl['halo_posts']
  for arr,name,dt in [(z,'height_mm.i32','<i4'),(shore,'shore_mm.i32','<i4'),(fixed,'fixed.u8','u1')]:
   a=np.fromfile(inp/d.name/'fields'/name,dt).reshape(nl['n'],nl['n'])[nh:nh+512,nh:nh+512];arr[yy:yy+512,xx:xx+512]=a/(1000 if name!='fixed.u8' else 1)
  for bit in [1,2,4,8,16,32,64,128,256]:flagcounts[str(bit)]=flagcounts.get(str(bit),0)+int(np.count_nonzero(fields['coast_flags'][core]&bit))
 for (x,y),(d,h) in wins.items():
  for dx,dy in [(1,0),(0,1)]:
   if (x+dx,y+dy) not in wins:continue
   e,hh=wins[x+dx,y+dy];a=h['lattice'];b=hh['lattice'];x0=max(a['gx0'],b['gx0']);y0=max(a['gy0'],b['gy0']);x1=min(a['gx0']+a['n'],b['gx0']+b['n']);y1=min(a['gy0']+a['n'],b['gy0']+b['n']);diff={}
   for f,g in zip(h['fields'],hh['fields']):
    aa=read(d,f,a['n'])[y0-a['gy0']:y1-a['gy0'],x0-a['gx0']:x1-a['gx0']];bb=read(e,g,b['n'])[y0-b['gy0']:y1-b['gy0'],x0-b['gx0']:x1-b['gx0']]
    diff[f['name']]=int(np.count_nonzero(aa.view('u1')!=bb.view('u1')))
   seams.append({'pair':[d.name,e.name],'overlap_posts':(x1-x0)*(y1-y0),'unequal_bytes':diff,'exact':not any(diff.values())})
 assert len(seams)==12 and all(s['exact'] for s in seams)
 counts=np.bincount(cls.ravel(),minlength=8);n=int(counts[1:].sum());shares={NAMES[i]:{'posts':int(counts[i]),'share_of_classified':float(counts[i]/n)} for i in range(1,8)}
 assert not np.any((cls>0)&(fixed!=0)), 'coast on water/channel post'
 assert not np.any(np.isin(cls,[5,6])&((z < -1.0)|(z>1.0))), 'flat or marsh outside tidal frame'
 assert not np.any(np.isin(cls,[5,6])&(sl>math.tan(math.radians(3))+1e-8)), 'flat or marsh exceeds local gradient'
 near=(fixed!=2)&(fixed!=4)&(shore>=0)&(shore<=32);ss=slope(z,32,32);valid=near&np.isfinite(ss)
 proxy={'land_shore_posts':int(valid.sum()),'steep_ge30_share':float(np.mean(ss[valid]>=math.tan(math.radians(30)))),'low_lt10_share':float(np.mean(ss[valid]<math.tan(math.radians(10))))}
 paint(z,cls,notes/f'{tag}_coast_hillshade.png',f'C7 seed 42 site {tag[1:]}: coast potential on native terrain',f'3 x 3 windows, 8.192 km native evaluation margin, {n:,} classified land posts. Seams 12/12 exact.',fixed=fixed)
 # Close crop selected mechanically around class activity, no selection of attractive outcome.
 cx,cy={'s1':(563,612),'s5':(668,273)}[tag];r=192 # Frozen arc1 crop for honest before/after
 paint(z[cy-r:cy+r,cx-r:cx+r],cls[cy-r:cy+r,cx-r:cx+r],notes/f'{tag}_coast_detail.png',f'C7 {tag}: 12.3 km coast detail','Same crop as arc1. Blue = native water. Tan = low native ground; includes constructed bars.',fixed=fixed[cy-r:cy+r,cx-r:cx+r])
 return {'seed':42,'input':str(inp),'output':str(out),'seams':seams,'shares':shares,'flags_half_open_core':flagcounts,'morphology_proxy':proxy,'fields':manifest}

def real_dem(root,notes):
 import rasterio
 from rasterio.merge import merge
 records={}
 for name,tags,bounds in [('kauai',['N21_00_W160','N22_00_W160'],(-159.85,21.8,-159.25,22.25)),('maui',['N20_00_W157','N20_00_W156'],(-156.75,20.5,-155.95,21.0))]:
  files=[next(root.glob('tiles/*'+tag+'*/*DEM.tif')) for tag in tags];src=[rasterio.open(p) for p in files]
  a,t=merge(src,bounds=bounds);z=a[0];lat=(bounds[1]+bounds[3])/2;dx=111195*math.cos(math.radians(lat))*t.a;dy=111195*abs(t.e)
  sea=z<=0;near=np.zeros_like(sea)
  for yy in [-1,0,1]:
   for xx in [-1,0,1]:near|=np.roll(np.roll(sea,yy,0),xx,1)
  ss=slope(z,dx,dy);valid=(z>0)&near&np.isfinite(ss);steep=ss>=math.tan(math.radians(30));low=ss<math.tan(math.radians(10));classes=np.zeros_like(z,np.uint8);classes[valid&steep]=1;classes[valid&low]=4;classes[valid&~steep&~low]=2
  records[name]={'files':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files],'bounds_lonlat':bounds,'spacing_m':[dx,dy],'coastal_land_posts':int(valid.sum()),'steep_ge30_share':float(np.mean(steep[valid])),'low_lt10_share':float(np.mean(low[valid])),'limitation':'GLO30 DSM and zero-elevation sea proxy, not observed beach/cliff ground truth. Low slope is only beach-compatible. Cropped island coverage; no grain or wave data.'}
  # raster north first, paint expects south first
  paint(np.flipud(z),np.flipud(classes),notes/f'{name}_dem_proxy.png',f'{name.title()}: GLO30 coastal slope proxy','Red >=30 deg; yellow <10 deg; purple intermediate. These are slope bins, not measured shore classes.',spacing=dx,dy=dy,real=True)
  for s in src:s.close()
 return records

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--outputs',type=P,required=True);ap.add_argument('--notes',type=P,required=True);ap.add_argument('--site1',type=P,required=True);ap.add_argument('--site5',type=P,required=True);ap.add_argument('--real-dem',type=P);a=ap.parse_args();a.notes.mkdir(parents=True,exist_ok=True)
 report={s:coast_site(inp,a.outputs/s,s,a.notes) for s,inp in [('s1',a.site1),('s5',a.site5)]}
 if a.real_dem:report['real_dem']=real_dem(a.real_dem,a.notes)
 (a.notes/'acceptance.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk in ['shares','morphology_proxy']} for k,v in report.items() if k!='real_dem'},indent=2))
if __name__=='__main__':main()
