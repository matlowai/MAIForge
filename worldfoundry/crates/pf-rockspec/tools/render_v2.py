#!/usr/bin/env python3
"""Side-by-side v1/v2 response and CPU-rendered material-coordinate closeups.
Consumes outputs only. PIL is the pre-existing local install; no new packages.
"""
import sys,runpy,json,math
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
p=Path(sys.argv[1]);root=Path(__file__).resolve().parent
# Retain the exact v1 diagnostic renderer and lighting as the sphere comparator.
old=runpy.run_path(str(root/'render_swatches.py'),run_name='__main__')
sphere=old['sphere'];font=old['font'];small=old['small'];large=old['large']
cat=json.loads((root.parent/'data/catalog_v2.json').read_text())
sets=[json.loads((p/f'swatches_v{v}.json').read_text())['samples'] for v in [1,2]]
lookup=[{(x['input']['material_ref']['spec_id'],x['input']['footprint_m'],x['input']['surface']['surface_id'].split('/')[1]):x['output'] for x in xs} for xs in sets]
weather={(x['rock'],x['scenario']):x['output'] for x in json.loads((p/'weathering_v2.json').read_text())['samples']}
sanity={x['rock']:x for x in json.loads((p/'sanity_v2.json').read_text())['samples']}
LW,CW,HH,RH=205,90,110,244;W=LW+18*CW;H=HH+len(cat['rocks'])*RH
sheet=Image.new('RGB',(W,H),(24,26,29));draw=ImageDraw.Draw(sheet)
draw.text((10,8),'pf-rockspec v2 | V1 preserved beside V2 | 43 conditional materials',font=large,fill='white')
draw.text((10,34),'Top rows: same sphere and light; evaluator footprints. Below: physical 0.1/2 m patches (128 pixels), then homogenized 50 m response.',font=font,fill='white')
draw.text((10,52),'Weathered top = explicit wet-temperate residue. Arid comparison below. These are optical diagnostics, not world-admitted rock or Unreal.',font=font,fill=(200,205,210))
for side in range(2):
 for j,f in enumerate([.1,2.,50.]):
  for k,h in enumerate(['dry','wet','rind']):
   x=LW+(side*9+j*3+k)*CW;draw.text((x+2,77),f'v{side+1} {f:g}m',font=font,fill='white');draw.text((x+2,93),h,font=small,fill=(210,210,210))
for i,r in enumerate(cat['rocks']):
 y=HH+i*RH;k=r['id'];s=sanity[k]
 if i%2==0:draw.rectangle((0,y,W,y+RH),fill=(30,32,35))
 draw.text((8,y+6),k,font=font,fill='white');draw.text((8,y+27),r['structure_v2']['kind'],font=small,fill=(200,205,210))
 draw.text((8,y+43),f"optical boundary eff {r['structure_v2']['optical_boundary_efficiency']['value']:.3g}",font=small,fill=(200,205,210))
 for side in range(2):
  for j,f in enumerate([.1,2.,50.]):
   for hidx,h in enumerate(['fresh_dry','fresh_wet','weathered_dry']):
    out=lookup[side][(k,f,h)];x=LW+(side*9+j*3+hidx)*CW
    sheet.paste(sphere(out).resize((88,48),Image.Resampling.LANCZOS),(x,y+3));lum=sum(a*b for a,b in zip(out['base_reflectance_linear'],[.2126,.7152,.0722]))
    draw.text((x+2,y+54),f'Y {lum:.3f}',font=small,fill='white');draw.text((x+2,y+68),f'r {out["roughness"]:.2f}',font=small,fill=(200,205,210))
 # Detail row spans BOTH halves, preserving all nine comparison columns above.
 for j,span in enumerate([.1,2.,50.]):
  x=LW+j*164;filename=p/f'meso_{k}_{span:g}.ppm'
  # Rust float display writes 0.1, 2, 50; labels explicitly distinguish filtering.
  image=Image.open(filename);sheet.paste(image,(x,y+105))
  draw.text((x,y+87),f'{span:g}m '+('cell average' if span==50 else 'physical patch'),font=small,fill='white')
 for j,scenario in enumerate(['wet_temperate','arid_varnish']):
  x=LW+520+j*150;out=weather[(k,scenario)];sheet.paste(sphere(out),(x,y+119));draw.text((x,y+96),scenario+' | 2m',font=small,fill='white')
  lum=sum(a*b for a,b in zip(out['base_reflectance_linear'],[.2126,.7152,.0722]));draw.text((x,y+202),f'Y {lum:.3f}',font=font,fill='white')
 x=LW+850;lines=[f"2m: v1 {s['v1_Y']:.3f} -> v2 {s['v2_Y_2m']:.3f}",f"Integrator band: {s['integrator_band']} | {s['integrator_comparison']}",f"Cited specimen spectral span: {s['cited_band']}",f"Source comparison: {s['literature_comparison']} (not calibration)","Full metric/source limitations: sanity_table_v2.md"]
 for j,line in enumerate(lines):draw.text((x,y+107+22*j),line,font=font,fill=(210,215,220))
 draw.line((0,y+RH-1,W,y+RH-1),fill=(70,73,77))
sheet.save(p/'contact_sheet_v2_beside_v1.png')
for page,start in enumerate(range(0,len(cat['rocks']),5),1):
 count=min(5,len(cat['rocks'])-start);img=Image.new('RGB',(W,HH+count*RH));img.paste(sheet.crop((0,0,W,HH)),(0,0));img.paste(sheet.crop((0,HH+start*RH,W,HH+(start+count)*RH)),(0,HH));img.save(p/f'comparison_page_{page:02}.png')
# A larger, easy-to-review plate of all 43 closeup rows, same source pixels.
MW,MH=640,55+43*166;meso=Image.new('RGB',(MW,MH),(24,26,29));d=ImageDraw.Draw(meso);d.text((10,10),'V2 structure: 0.1m | 2m | 50m cell average',font=large,fill='white')
for i,r in enumerate(cat['rocks']):
 y=55+i*166;d.text((10,y),r['id'],font=font,fill='white')
 for j,span in enumerate([.1,2.,50.]):meso.paste(Image.open(p/f'meso_{r["id"]}_{span:g}.ppm'),(180+j*148,y+23))
meso.save(p/'mesostructure_all_v2.png')
for page,start in enumerate(range(0,43,15),1):
 count=min(15,43-start);im=Image.new('RGB',(MW,55+count*166));im.paste(meso.crop((0,0,MW,55)),(0,0));im.paste(meso.crop((0,55+start*166,MW,55+(start+count)*166)),(0,55));im.save(p/f'mesostructure_page_{page}.png')
(p/'render_manifest_v2.json').write_text(json.dumps(dict(width=W,height=H,rocks=43,comparison_swatches=774,closeups=129,exposure_swatches=86,pages=9,source='CPU evaluator output only; no authored textures or RGB rock presets',filter='0.1/2m image pixel footprint=span/128; 50m image is uniform 50m average; comparison sphere footprint=column label',renderer='v1 sphere lighting retained; texture plane uses same GGX family and bounded optical normals',limits='optical normals only, no displacement, terrain joints, silhouettes, collision, or backend parity'),indent=2)+'\n')
print('43 comparison rows + 43 closeup rows, 9 pages and 3 texture pages')
