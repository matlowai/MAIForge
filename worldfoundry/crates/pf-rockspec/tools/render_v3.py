#!/usr/bin/env python3
"""Render a labelled contact sheet from CPU v3 output only. Existing Pillow, offline."""
import json,sys,math
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
p=Path(sys.argv[1]);rows=json.loads((p/'samples_v3.json').read_text())['samples'];lookup={r['id']:r for r in rows}
fontpath='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
font=ImageFont.truetype(fontpath,14);small=ImageFont.truetype(fontpath,12);large=ImageFont.truetype(fontpath,21)
def lum(rgb):return sum(a*b for a,b in zip(rgb,[.2126,.7152,.0722]))
def srgb(v):v=max(0,min(1,v));return round(255*(12.92*v if v<=.0031308 else 1.055*v**(1/2.4)-.055))
def sphere(out,n=112):
 im=Image.new('RGB',(n,n),(24,27,31));pix=im.load();l=[-.3,-.4,math.sqrt(.75)]
 for y in range(n):
  for x in range(n):
   a=(x+.5-n/2)/(n*.46);b=(y+.5-n/2)/(n*.46);rr=a*a+b*b
   if rr>=1:continue
   z=math.sqrt(1-rr);nl=max(0,a*l[0]+b*l[1]+z*l[2])
   # One fixed diffuse sphere. No invented rock-specific tint or fitted gloss.
   pix[x,y]=tuple(srgb(v*(.15+.85*nl)) for v in out['diffuse_energy_linear'])
 return im
panels=[['granite.i_type_pink','granite.i_type_white','granite.s_type_pink','granite.s_type_white'],['chert.flint','chert.white','ignimbrite.unwelded','ignimbrite.welded'],['andesite.arc_hornblende','andesite.two_pyroxene','dacite.arc_hornblende','dacite.two_pyroxene'],['basalt_tholeiitic.lip','basalt_alkali.island_cap','shale.grey','shale.organic_rich']]
CW,RH,HH=320,235,100;im=Image.new('RGB',(4*CW,HH+4*RH),(24,27,31));d=ImageDraw.Draw(im)
d.text((12,8),'pf-rockspec v3 | Body realizations side by side',font=large,fill='white')
d.text((12,38),'Same light: 2 m filtered diffuse sphere (left), CPU material patch spanning 0.1 m (right).',font=font,fill='white')
d.text((12,60),'Conditional offline priors. V2 optical mixing retained. Not Unreal, calibrated spectra, or an admitted world.',font=font,fill=(235,195,125))
for row,ids in enumerate(panels):
 for col,id in enumerate(ids):
  r=lookup[id];x=col*CW;y=HH+row*RH;out=r['mean_2m'];d.rectangle((x,y,x+CW-1,y+RH-1),outline=(70,74,80));d.text((x+9,y+7),id,font=font,fill='white');d.text((x+9,y+28),r['setting'],font=small,fill=(190,200,210));im.paste(sphere(out),(x+9,y+50));im.paste(Image.open(p/r['file']).resize((112,112),Image.Resampling.NEAREST),(x+153,y+50));d.text((x+10,y+167),f"Y={lum(out['base_reflectance_linear']):.4f}  F0={sum(out['f0_linear'])/3:.3f}",font=font,fill='white');d.text((x+10,y+188),f"rho={out['dry_bulk_density_kg_m3']:.0f} kg/m3",font=small,fill=(205,210,220));params=r['parameters'];label=''
  if params.get('welding_grade') is not None:label=f"welding={params['welding_grade']:.1f}; fused-glass transport"
  elif params.get('organic_mass_fraction') is not None:label=f"organic dry mass={params['organic_mass_fraction']:.2f}"
  elif params.get('plagioclase_an_fraction') is not None:label=f"plagioclase An={params['plagioclase_an_fraction']:.2f}"
  d.text((x+10,y+208),label,font=small,fill=(205,210,220))
im.save(p/'contact_sheet_v3.png')
# Every facies reference plus all added realizations as a browsable diagnostic index.
W=4*CW;H=HH+math.ceil(len(rows)/4)*RH;full=Image.new('RGB',(W,H),(24,27,31));fd=ImageDraw.Draw(full);fd.text((12,12),'V3 full catalogue | 67 conditional realizations',font=large,fill='white')
for i,r in enumerate(rows):
 x=(i%4)*CW;y=HH+(i//4)*RH;out=r['mean_2m'];fd.text((x+8,y+6),r['id'],font=font,fill='white');full.paste(sphere(out),(x+8,y+45));full.paste(Image.open(p/r['file']).resize((112,112)),(x+150,y+45));fd.text((x+8,y+170),f"Y={lum(out['base_reflectance_linear']):.4f}; {out['dry_bulk_density_kg_m3']:.0f} kg/m3",font=font,fill='white')
full.save(p/'all_realizations_v3.png')
(p/'render_manifest_v3.json').write_text(json.dumps(dict(panels=16,all_realizations=len(rows),source='samples_v3.json and CPU ppm patches',patch_span_m=.1,pixels_per_patch=96,lighting='fixed diffuse directional+ambient; no exposure fitting',optical_mixing='unchanged v2; F3 not implemented'),indent=2)+'\n')
print(p/'contact_sheet_v3.png')
