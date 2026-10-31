#!/usr/bin/env python3
"""Lossless PNG delivery of the CPU 1024px 2m panels, with navigable index."""
import argparse,json
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
ap=argparse.ArgumentParser();ap.add_argument('detail',type=Path);ap.add_argument('run',type=Path);a=ap.parse_args()
metrics=json.loads((a.detail/'mesostructure_v2.json').read_text())
lines=['# Full-resolution 2 m structure panels','','CPU material-plane optical diagnostics, 1024 x 1024; physical pixel footprint 1.953125 mm. No image enlargement or generated geology. Fine structure below this support is filtered. These are not mesh displacement or natural-rock photographs.','']
for row in metrics['samples']:
 if row['span_m']!=2:continue
 f=a.detail/row['file'];im=Image.open(f);assert im.size==(1024,1024)
 dest=f.with_suffix('.png');im.save(dest)
 lines.append(f'- [{row["rock"]}]({dest.resolve()}): Y mean {row["Y_mean"]:.5f}, spatial standard deviation {row["Y_std"]:.5f}.')
(a.run/'DETAIL_INDEX.md').write_text('\n'.join(lines)+'\n')
# Review crops at native source pixels, paired with the small 0.1m patch.
keys=['andesite','dacite','scoria','pumice','granite','conglomerate','siltstone','shale','schist','gneiss','gravel','till']
font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',15)
plate=Image.new('RGB',(800,50+len(keys)*276),(24,26,29));d=ImageDraw.Draw(plate);d.text((10,10),'0.1m patch (2x display) | native-pixel crop from 2m panel',font=font,fill='white')
for i,k in enumerate(keys):
 y=50+i*276;d.text((8,y),k,font=font,fill='white')
 im=Image.open(a.detail/f'meso_{k}_0.1.ppm');plate.paste(im.resize((256,256),Image.Resampling.NEAREST),(180,y+18))
 im=Image.open(a.detail/f'meso_{k}_2.ppm');plate.paste(im.crop((384,384,640,640)),(470,y+18))
plate.save(a.run/'evidence/detail_review_v2.png')
for page,start in enumerate(range(0,len(keys),4),1):
 im=Image.new('RGB',(800,50+4*276));im.paste(plate.crop((0,0,800,50)),(0,0));im.paste(plate.crop((0,50+start*276,800,50+(start+4)*276)),(0,50));im.save(a.run/f'evidence/detail_review_page_{page}.png')
print('43 native-resolution PNG panels; 12 requested motif review rows')
