#!/usr/bin/env python3
"""Render ONLY exported CPU outputs, not a second material evaluator.
Neutral GGX spheres reveal reflectance and lobe width. No AI images, textures,
key-based colours, geology generation, downloads or external assets.
"""
import argparse,json,math
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
ap=argparse.ArgumentParser();ap.add_argument('evidence',type=Path);args=ap.parse_args()
root=Path(__file__).resolve().parents[1]
catalog=json.loads((root/'data/catalog_v1.json').read_text())
data=json.loads((args.evidence/'swatches_v1.json').read_text())
lookup={x['input']['id']:x for x in data['samples']}
font_path='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
font=ImageFont.truetype(font_path,12);small=ImageFont.truetype(font_path,10);large=ImageFont.truetype(font_path,20)
CW,CH,LW,HH=155,102,252,112
W=LW+9*CW;H=HH+len(catalog['rocks'])*CH
sheet=Image.new('RGB',(W,H),(24,26,29));draw=ImageDraw.Draw(sheet)
draw.text((12,10),'pf-rockspec v1 | 43 rock priors | CPU reflectance / GGX reference swatches',font=large,fill='white')
draw.text((12,38),'Fixed neutral light, identical smooth geometry, sRGB display. Y labels are LINEAR base-reflectance luminance, not displayed brightness.',font=font,fill=(210,210,210))
draw.text((12,57),'Columns give evaluator filter footprint: grains below it are averaged. Weathered states are explicit diagnostic histories, not predicted exposure ages.',font=font,fill=(210,210,210))
for j,f in enumerate([.1,2.,50.]):
 for k,h in enumerate(['Fresh dry','Fresh wet','Weathered dry']):
  x=LW+(j*3+k)*CW; draw.text((x+4,79),f'{f:g} m | {h}',font=font,fill='white')
L=(-.35,-.45,math.sqrt(1-.35**2-.45**2));hv=(L[0],L[1],L[2]+1);hn=math.sqrt(sum(x*x for x in hv));HV=tuple(x/hn for x in hv)
clip_count=0;pixel_count=0
# Geometry is shared. Deliberately no invented mineral speckle or tint detail.
geom=[]
for y in range(73):
 for x in range(143):
  nx=(x-71)/59;ny=(y-34)/32;rr=nx*nx+ny*ny
  if rr>=1:continue
  nz=math.sqrt(1-rr);nl=max(0.,nx*L[0]+ny*L[1]+nz*L[2]);nh=max(0.,nx*HV[0]+ny*HV[1]+nz*HV[2]);vh=HV[2]
  geom.append((x,y,nz,nl,nh,vh))
def srgb(x):
 x=max(0,min(1,x));return int(round(255*(12.92*x if x<=.0031308 else 1.055*x**(1/2.4)-.055)))
def sphere(out):
 global clip_count,pixel_count
 img=Image.new('RGB',(143,73),(32,34,37));pix=img.load()
 alpha=max(.002,out['roughness']**2);a2=alpha*alpha
 for x,y,nv,nl,nh,vh in geom:
  d=a2/(math.pi*(nh*nh*(a2-1)+1)**2)
  def g1(n): return 2*n/(n+math.sqrt(a2+(1-a2)*n*n)) if n>0 else 0.
  g=g1(nv)*g1(nl)
  rgb=[]
  for ch in range(3):
   f0=out['f0_linear'][ch];f=f0+(1-f0)*(1-vh)**5
   # Unit-white irradiance convention: E_direct=pi; neutral low ambient.
   value=out['base_reflectance_linear'][ch]*(1-f)*nl + math.pi*d*g*f/(4*max(nv,1e-9)) + .08*out['diffuse_energy_linear'][ch]
   if value>1:clip_count+=1
   pixel_count+=1;rgb.append(srgb(value))
  pix[x,y]=tuple(rgb)
 return img
for row,r in enumerate(catalog['rocks']):
 y=HH+row*CH
 if row%2==0:draw.rectangle((0,y,W,y+CH),fill=(29,31,34))
 draw.text((10,y+12),r['id'],font=font,fill='white')
 mode=sorted(r['mode'],key=lambda m:-m['mass_fraction']['value'])[:3]
 for i,m in enumerate(mode):
  draw.text((10,y+34+15*i),f"{m['mineral']} {m['mass_fraction']['value']:.0%} mass",font=small,fill=(188,194,200))
 for j,f in enumerate([.1,2.,50.]):
  for k,h in enumerate(['fresh_dry','fresh_wet','weathered_dry']):
   cell=lookup[f'{r["id"]}/{h}/{f}'];out=cell['output'];x=LW+(j*3+k)*CW
   sheet.paste(sphere(out),(x+5,y+3))
   rgb=out['base_reflectance_linear'];lum=sum(a*b for a,b in zip(rgb,[.2126,.7152,.0722]))
   draw.text((x+7,y+80),f'Y {lum:.3f}   r {out["roughness"]:.2f}',font=font,fill=(230,230,230))
   if out['coating'] and out['coating']['optical_weight']>0:draw.rectangle((x+143,y+80,x+147,y+90),fill=(180,130,70))
sheet.save(args.evidence/'contact_sheet_v1.png')
# Readable pages preserve source pixels; crops are layout outputs, not edits to an image.
for page,start in enumerate(range(0,len(catalog['rocks']),15),1):
 count=min(15,len(catalog['rocks'])-start);img=Image.new('RGB',(W,HH+CH*count),(24,26,29));img.paste(sheet.crop((0,0,W,HH)),(0,0));img.paste(sheet.crop((0,HH+start*CH,W,HH+(start+count)*CH)),(0,HH));img.save(args.evidence/f'contact_sheet_page_{page}.png')
meta=dict(width=W,height=H,swatches=387,renderer='fixed smooth ellipsoid, neutral direct GGX + 0.08 diffuse ambient, no key textures',light=L,colour='linear to IEC sRGB, display clamp at 1; Y is unclamped prelight base reflectance',clipped_channel_values=clip_count,evaluated_channel_values=pixel_count,footprint_semantics='CPU mean filter footprint, not pixel footprint or a textured physical rock photograph',sparkle='CPU energy fraction exported; no unfiltered sparkle dots drawn at averaged footprints',normal_detail='unresolved variance enters CPU roughness; resolved microgeometry not rendered on smooth diagnostic shape')
(args.evidence/'render_manifest_v1.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta))
