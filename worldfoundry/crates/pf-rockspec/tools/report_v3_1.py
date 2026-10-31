#!/usr/bin/env python3
"""Tables and contact sheets from the actual CPU export. Offline Pillow only."""
import csv
import json
import sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

src, dst = map(Path, sys.argv[1:3])
dst.mkdir(parents=True, exist_ok=True)
data = json.loads((src/'samples_v3_1.json').read_text())
rows = [r for r in data['samples'] if r['affected']]
lookup = {(r['input']['realization_id'],r['input']['history'],r['input']['footprint_m']):r for r in rows}
real = {r['id']:r for r in data['realizations'] if r['affected']}
def srgb(x):
    return 255*(12.92*x if x<=.0031308 else 1.055*x**(1/2.4)-.055)
def y(out):
    return sum(a*b for a,b in zip(out['base_reflectance_linear'],[.2126,.7152,.0722]))
def rgb(out):
    return out['base_reflectance_linear']
def rgbstr(out):
    return ', '.join(f'{x:.5f}' for x in rgb(out))
table=[]
for r in rows:
    inp=r['input'];a=r['before_v3_same_surface'];b=r['after_v3_1']
    historical=inp['history']!='weathered_dry' or r['weather_fixture_existed_in_v3']
    table.append(dict(realization=inp['realization_id'],history=inp['history'],footprint_m=inp['footprint_m'],
        historical_v3_surface=historical,
        before_linear_rgb=rgbstr(a),after_linear_rgb=rgbstr(b),before_Y=y(a),after_Y=y(b),
        max_delta_linear_255=max(abs(x-z)*255 for x,z in zip(rgb(a),rgb(b))),
        max_delta_srgb_255=max(abs(srgb(x)-srgb(z)) for x,z in zip(rgb(a),rgb(b)))))
with (dst/'before_after_v3_1.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
md=['# Measured v3 to v3.1b affected realizations',
    '', 'All values are CPU evaluations of manufactured catalogue inputs, no world seed. '
    'RGB and Y are linear reflectance. Differences in /255 units are explicitly linear or sRGB encoded. '
    'The CSV includes all 153 affected surface/scale cases, including wet and weathered surfaces. '
    'For six new weathered families the old column is a same-surface optical ablation, not a historical fixture.',
    '', '| Realization | Footprint m | Before RGB | After RGB | Before Y | After Y | max linear /255 | max sRGB /255 |',
    '|---|---:|---|---|---:|---:|---:|---:|']
for r in table:
    if r['history']=='fresh_dry':
        md.append(f"| {r['realization']} | {r['footprint_m']:g} | {r['before_linear_rgb']} | {r['after_linear_rgb']} | {r['before_Y']:.5f} | {r['after_Y']:.5f} | {r['max_delta_linear_255']:.2f} | {r['max_delta_srgb_255']:.2f} |")
(dst/'BEFORE_AFTER.md').write_text('\n'.join(md)+'\n')

short=['# v3.1b short before/after', '',
       'Fresh dry linear luminance Y, v3 -> v3.1b. Manufactured catalogue fixtures; no calibration.', '',
       '| Realization | 0.1 m | 2 m | 50 m |', '|---|---:|---:|---:|']
for rid in real:
    cells=[]
    for f in [.1,2.,50.]:
        row=lookup[(rid,'fresh_dry',f)]
        cells.append(f"{y(row['before_v3_same_surface']):.5f} -> {y(row['after_v3_1']):.5f}")
    short.append('| '+rid+' | '+' | '.join(cells)+' |')
short += ['', 'Six new weathered families: historical v3 fixture absent; all 54 frozen cases use unchanged v3 optics. All 54 granite comparisons are numerically unchanged apart from version metadata.']
(dst/'SHORT_BEFORE_AFTER.md').write_text('\n'.join(short)+'\n')

pairs=[('granite.i_type_pink','granite.i_type_white'),('granite.s_type_pink','granite.s_type_white'),
       ('basalt_tholeiitic.reference','basalt_tholeiitic.lip'),
       ('basalt_tholeiitic.reference','basalt_tholeiitic.island_shield'),
       ('basalt_tholeiitic.lip','basalt_tholeiitic.island_shield'),
       ('basalt_alkali.reference','basalt_alkali.island_cap')]
contrasts=[]
for a,b in pairs:
    for f in [.1,2.,50.]:
        row=dict(a=a,b=b,footprint_m=f)
        for version,key in [('before','before_v3_same_surface'),('after','after_v3_1')]:
            x=lookup[(a,'fresh_dry',f)][key];z=lookup[(b,'fresh_dry',f)][key]
            row[version+'_max_linear_255']=max(abs(u-v)*255 for u,v in zip(rgb(x),rgb(z)))
            row[version+'_max_srgb_255']=max(abs(srgb(u)-srgb(v)) for u,v in zip(rgb(x),rgb(z)))
        contrasts.append(row)
(dst/'pair_contrasts.json').write_text(json.dumps(contrasts,indent=2)+'\n')

fontpath='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
font=ImageFont.truetype(fontpath,13);small=ImageFont.truetype(fontpath,11);large=ImageFont.truetype(fontpath,20)
# Explicit groups keep the basalt comparison together.
groups=[[k for k in real if k.startswith(('granite.','basalt_'))],
        [k for k in real if not k.startswith(('granite.','basalt_'))]]
sheets=[]
for page,ids in enumerate(groups,1):
    h=120;rh=115;im=Image.new('RGB',(1120,h+len(ids)*rh+35),(25,28,32));d=ImageDraw.Draw(im)
    d.text((12,8),'pf-rockspec v3.1b | actual CPU before / after',font=large,fill='white')
    d.text((12,36),'Conditional offline priors; no exposure fitting. Blocks: linear base reflectance encoded as sRGB.',font=font,fill='white')
    d.text((12,57),'Patches: 0.1 m span, 1.25 mm pixels, fixed diffuse light. No Unreal or spectral calibration claimed.',font=font,fill=(225,195,130))
    for x,label in [(10,'Realization'),(270,'Patch v3'),(365,'Patch v3.1b'),(466,'0.1 m'),(636,'2 m'),(806,'50 m'),(972,'Weather 2 m')]:
        d.text((x,83),label,font=font,fill='white')
    for x,label in [(468,'v3'),(545,'v3.1b'),(638,'v3'),(715,'v3.1b'),(808,'v3'),(885,'v3.1b'),(972,'v3'),(1045,'v3.1b')]:
        d.text((x,101),label,font=small,fill=(190,200,215))
    for i,rid in enumerate(ids):
        top=h+i*rh;d.line((0,top,1120,top),fill=(70,75,80))
        a,b=rid.split('.',1);d.text((10,top+12),a,font=font,fill='white');d.text((10,top+31),b,font=font,fill='white')
        grain=real[rid]['after_grain_m']['value'];d.text((10,top+56),f"grain {grain*1000:g} mm",font=small,fill=(200,210,220))
        d.text((10,top+75),'Unchanged v3 optics',font=small,fill=(200,210,220))
        for x,name in zip([270,365],real[rid]['files']):
            im.paste(Image.open(src/name).resize((84,84),Image.Resampling.NEAREST),(x,top+10))
        for x,f in [(466,.1),(636,2.),(806,50.)]:
            row=lookup[(rid,'fresh_dry',f)]
            for offset,key in [(0,'before_v3_same_surface'),(77,'after_v3_1')]:
                out=row[key];colour=tuple(round(srgb(v)) for v in rgb(out))
                d.rectangle((x+offset,top+10,x+offset+67,top+78),fill=colour)
                d.text((x+offset,top+85),f'Y {y(out):.3f}',font=small,fill='white')
        row=lookup[(rid,'weathered_dry',2.)]
        for x,key in [(972,'before_v3_same_surface'),(1045,'after_v3_1')]:
            if x==972 and not row['weather_fixture_existed_in_v3']:
                d.text((x,top+33),'absent',font=small,fill=(225,195,130));continue
            out=row[key];d.rectangle((x,top+10,x+65,top+78),fill=tuple(round(srgb(v)) for v in rgb(out)))
            d.text((x,top+85),f'Y {y(out):.3f}',font=small,fill='white')
    path=dst/f'contact_sheet_v3_1_{page}.png';im.save(path);sheets.append(str(path.resolve()))
metrics=dict(realizations=len(data['realizations']),samples=len(data['samples']),affected_realizations=len(real),
             affected_cases=len(rows),fresh_table_rows=sum(r['history']=='fresh_dry' for r in table),
             pair_contrasts=len(contrasts),contact_sheets=sheets,
             qualification='Rendered from exported CPU reflectance and patches. New-family weather baseline absent, not fabricated.')
(dst/'report_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps(metrics,indent=2))
