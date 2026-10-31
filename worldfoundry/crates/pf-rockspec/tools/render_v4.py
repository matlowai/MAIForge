#!/usr/bin/env python3
"""Offline CPU diagnostic contact sheets. Reads actual evaluator outputs only."""
import json, math, sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

p = Path(sys.argv[1])
rows = json.loads((p/'samples_v4.json').read_text())['samples']
fontpath = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
font = ImageFont.truetype(fontpath, 14)
small = ImageFont.truetype(fontpath, 12)
large = ImageFont.truetype(fontpath, 22)
clipped = 0
channels = 0
def y(out): return sum(a*b for a,b in zip(out['base_reflectance_linear'], [.2126,.7152,.0722]))
def srgb(x):
    global clipped, channels
    channels += 1
    clipped += x > 1
    x = max(0,min(1,x))
    return round(255*(12.92*x if x <= .0031308 else 1.055*x**(1/2.4)-.055))
def sphere(out, n=96):
    im=Image.new('RGB',(n,n),(24,27,31)); pix=im.load()
    light=[-.35,-.45,.82]; mag=math.sqrt(sum(v*v for v in light)); light=[v/mag for v in light]
    h=[light[0],light[1],light[2]+1];mag=math.sqrt(sum(v*v for v in h));h=[v/mag for v in h]
    alpha2=max(out['roughness'],.045)**4
    for py in range(n):
        for px in range(n):
            a=(px+.5-n/2)/(n*.46); b=(py+.5-n/2)/(n*.46)
            if a*a+b*b>=1: continue
            normal=[a,b,math.sqrt(1-a*a-b*b)]; nv=max(normal[2],.001)
            nl=max(0,sum(a*b for a,b in zip(normal,light))); nh=max(0,sum(a*b for a,b in zip(normal,h)))
            d=alpha2/(math.pi*(nh*nh*(alpha2-1)+1)**2)
            smith=lambda n:2*n/(n+math.sqrt(alpha2+(1-alpha2)*n*n))
            spec=math.pi*d*smith(nv)*smith(nl)/(4*nv)
            rgb=[]
            for base,f0 in zip(out['diffuse_energy_linear'],out['f0_linear']):
                fresnel=f0+(1-f0)*(1-h[2])**5
                rgb.append(srgb(base*(.15+.85*nl)+.85*spec*fresnel))
            pix[px,py]=tuple(rgb)
    return im

W,RH,HH=1200,132,112
pages=[]
for page in range(math.ceil(len(rows)/8)):
    group=rows[page*8:(page+1)*8]
    im=Image.new('RGB',(W,HH+len(group)*RH),(24,27,31));d=ImageDraw.Draw(im)
    d.text((12,8),f'pf-rockspec v3 / v4 | realization comparison {page+1}',font=large,fill='white')
    d.text((12,38),'Actual CPU material outputs. Same light, sRGB encoding, 2 m spheres and 0.1 m patches.',font=font,fill='white')
    d.text((12,60),'Offline priors, no spectral calibration or Unreal acceptance. Bright coarse-rock overshoot remains visible.',font=font,fill=(240,195,125))
    for x,label in [(12,'Realization'),(310,'v3 sphere / patch'),(534,'v4 sphere / patch'),(764,'Linear Y; shine and morphology')]:
        d.text((x,88),label,font=font,fill=(190,200,215))
    for i,r in enumerate(group):
        yy=HH+i*RH;d.line((0,yy,W,yy),fill=(66,72,80))
        s=r['scales'][1];d.text((12,yy+12),r['id'],font=font,fill='white')
        d.text((12,yy+40),f"grain {r['mixing']['physical_grain_m']*1e3:.4g} mm",font=small,fill=(190,200,215))
        for version,x in [(3,310),(4,534)]:
            out=s[f'v{version}'];im.paste(sphere(out),(x,yy+10))
            im.paste(Image.open(p/f"{r['id']}.v{version}.ppm"),(x+106,yy+10))
        d.text((764,yy+12),f"Y  {y(s['v3']):.4f} -> {y(s['v4']):.4f}",font=font,fill='white')
        d.text((764,yy+36),f"F0 {sum(s['v3']['f0_linear'])/3:.3f} -> {sum(s['v4']['f0_linear'])/3:.3f}",font=small,fill='white')
        d.text((764,yy+57),f"roughness {s['v3']['roughness']:.3f} -> {s['v4']['roughness']:.3f}",font=small,fill='white')
        if r['morphology']:
            m=r['morphology'];d.text((764,yy+80),f"Coat {m['thickness_m']*1e6:.2f} um / {m['coverage']:.0%} grains",font=small,fill=(240,195,125))
            d.text((764,yy+100),f"Dispersed-coat ablation Y {y(s['v4_without_coatings']):.4f}",font=small,fill=(240,195,125))
        else:
            d.text((764,yy+80),'Areal weights '+', '.join(f'{v:.2f}' for v in r['mixing']['areal_weight']),font=small,fill=(190,200,215))
    name=f'comparison_v4_page_{page+1:02}.png'; im.save(p/name);pages.append(name)
# One full sheet preserves the requested realization rows; pages are easier to view.
full=Image.new('RGB',(W,HH+RH*len(rows)),(24,27,31))
for i,name in enumerate(pages):
    im=Image.open(p/name)
    if i==0:full.paste(im,(0,0))
    else:full.paste(im.crop((0,HH,W,im.height)),(0,HH+i*8*RH))
full.save(p/'contact_sheet_v4_vs_v3.png')
(p/'render_manifest_v4.json').write_text(json.dumps(dict(realizations=len(rows),patches=len(rows)*2,pages=pages,
    full_sheet='contact_sheet_v4_vs_v3.png',shader='shared diffuse energy + fixed-light GGX sphere; diffuse structure patches',
    display_clipped_sphere_channels=clipped,sphere_channels=channels,exposure=1,
    qualification='CPU diagnostics, same admitted composition refs. Display highlight clipping is recorded, not material energy gain.'),indent=2)+'\n')
print(f'{len(rows)} realization rows; {len(pages)} pages; clipped sphere channels {clipped}/{channels}')
