#!/usr/bin/env python3
"""Label real CPU evaluator patches with existing Pillow. Offline, no fitting."""
import json
import sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

source = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else source
out.mkdir(parents=True, exist_ok=True)
data = json.loads((source / 'samples_v3_2b.json').read_text())
rows = {r['id']: r for r in data['rows']}
font_path = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
font = ImageFont.truetype(font_path, 12)
small = ImageFont.truetype(font_path, 10)
large = ImageFont.truetype(font_path, 18)
scales = [(0.001, '1 mm'), (0.1, '0.1 m'), (2, '2 m'), (10, '10 m')]
selected = ['schist.reference', 'phyllite.reference', 'gneiss.reference',
            'banded_iron_formation.reference', 'ignimbrite.welded',
            'ignimbrite.unwelded', 'pumice.reference', 'scoria.reference',
            'granite.reference', 'shale.reference', 'sandstone_quartz.reference',
            'basalt_tholeiitic.reference']

def patch(row, scale, version):
    return next(p for p in row['patches'] if p['span_m'] == scale and p['version'] == version)

def make_sheet(ids, channel, name):
    label, size, gap, top, height = 196, 128, 12, 112, 158
    width = label + 4 * (size * 2 + gap)
    im = Image.new('RGB', (width, top + len(ids) * height + 28), (24, 27, 31))
    d = ImageDraw.Draw(im)
    d.text((12, 10), 'pf-rockspec v3.2b structure | same body, position and light', font=large, fill='white')
    d.text((12, 36), 'Each scale: v3.2 left | opt-in v3.2b right. Patch width is physical size; pixel footprint = width / 128.', font=font, fill='white')
    d.text((12, 55), 'Manufactured catalogue fixtures. CPU shared evaluator. No UE parity, geological admission or calibrated spectra claimed.', font=font, fill=(235, 195, 125))
    d.text((12, 75), f'Channel: {channel}. Fixed exposure. Optical traces are not displacement, fractures or mineral stock.', font=font, fill=(195, 205, 215))
    key = {'GGX + diffuse': 'file', 'base reflectance (sRGB display)': 'albedo_file', 'material normal': 'normal_file'}[channel]
    for k, (_, scale_label) in enumerate(scales):
        x = label + k * (2 * size + gap)
        d.text((x, 93), f'{scale_label}     v3.2                       v3.2b', font=font, fill='white')
    for i, ident in enumerate(ids):
        y = top + i * height
        row = rows[ident]
        d.text((8, y + 6), ident.replace('.reference', ''), font=font, fill='white')
        rec = row['portable_structure']['record']
        fabric = rec['fabric']
        d.text((8, y + 27), fabric['kind'], font=small, fill=(195, 205, 215))
        spacing = fabric['resolved_spacing_m']
        if spacing:
            d.text((8, y + 45), f"fabric spacing {spacing['value']:g} m", font=small, fill=(195, 205, 215))
        if rec['foliated_cells']:
            d.text((8, y + 63), f"cell aspect {rec['fabric_cell_aspect']['value']:g}:1", font=small, fill=(195, 205, 215))
        if rec['welded_fiamme']:
            d.text((8, y + 63), f"lens flattening {rec['fiamme_flattening']['value']:g}:1", font=small, fill=(195, 205, 215))
        for k, (scale, _) in enumerate(scales):
            for v, version in enumerate(['v3_2', 'v3_2b']):
                x = label + k * (2 * size + gap) + v * size
                p = patch(row, scale, version)
                im.paste(Image.open(source / p[key]), (x, y))
                lum = sum(a * b for a, b in zip(p['mean_reflectance'], [.2126, .7152, .0722]))
                d.text((x + 2, y + size + 2), f'mean Y {lum:.4f}', font=small, fill=(195, 205, 215))
        d.line((0, y + height - 6, width, y + height - 6), fill=(60, 65, 70))
    d.text((10, im.height - 22), 'Finite patches need not have the bulk mean. Exact far-footprint equality and statistical mean checks are separate gates.', font=font, fill='white')
    im.save(out / name)

selected = [i for i in selected if i in rows]
make_sheet(selected, 'GGX + diffuse', 'contact_sheet_v3_2_vs_v3_2b.png')
make_sheet(selected, 'base reflectance (sRGB display)', 'contact_sheet_albedo_v3_2_vs_v3_2b.png')
make_sheet(selected, 'material normal', 'contact_sheet_normals_v3_2_vs_v3_2b.png')
for start in range(0, len(rows), 12):
    make_sheet(list(rows)[start:start + 12], 'GGX + diffuse', f'all_realizations_{start // 12 + 1:02}.png')

mica_ids = ['schist.reference', 'phyllite.reference', 'granite.reference']
im = Image.new('RGB', (1100, 110 + len(mica_ids) * 300), (24, 27, 31))
d = ImageDraw.Draw(im)
d.text((12, 10), 'Mica: one normal per pixel | foliation-parallel faces', font=large, fill='white')
d.text((12, 38), 'v3.2 and v3.2b use the same single-normal GGX. V3.2b selects bounded tilted facets by material cell.', font=font, fill='white')
d.text((12, 60), 'No additive sparkle energy. White mask shows selected facets; coloured image shows the actual material normal.', font=font, fill=(235, 195, 125))
for i, ident in enumerate(mica_ids):
    y = 110 + i * 300
    row = rows[ident]
    old, new = row['mica']
    d.text((8, y - 22), f"{ident} | width {rows[ident]['mica'][0]['span_m'] * 1000:g} mm", font=font, fill='white')
    for k, (p, field, title) in enumerate([(old, 'file', 'v3.2'), (new, 'file', 'v3.2b'), (new, 'normal_file', 'v3.2b normal'), (new, 'mask_file', 'v3.2b facet mask')]):
        x = 8 + k * 272
        im.paste(Image.open(source / p[field]), (x, y))
        d.text((x, y + 259), title, font=font, fill='white')
    d.text((820, y + 275), f"selected {new['flake_pixels']}/65536 pixels", font=small, fill='white')
im.save(out / 'mica_one_normal_v3_2_vs_v3_2b.png')
(out / 'render_manifest.json').write_text(json.dumps({
    'input': str(source.resolve()), 'rows': len(rows), 'selected_rows': selected,
    'scales_m': [x[0] for x in scales], 'n': 128,
    'channel_sheets': 3, 'atlas_pages': (len(rows) + 11) // 12, 'mica_sheets': 1,
    'lighting': 'fixed GGX + diffuse; no per-material fit; same one-normal path',
    'scope': 'manufactured shared-crate diagnostic, not an engine screenshot'
}, indent=2) + '\n')
print(out / 'contact_sheet_v3_2_vs_v3_2b.png')

# An explicitly labelled magnification of existing samples, not a finer-footprint
# render or exposure adjustment. It makes the thin fiamme profile inspectable.
r = rows['ignimbrite.welded']
detail = Image.new('RGB', (1584, 600), (24, 27, 31))
dd = ImageDraw.Draw(detail)
dd.text((12, 10), 'Welded ignimbrite | 0.1 m patch | existing 128 pixels enlarged 4x', font=large, fill='white')
dd.text((12, 38), 'Nearest-neighbour display only. Same 0.78125 mm sample footprint, fixed exposure and lighting.', font=font, fill='white')
for k, (version, field, title) in enumerate([
    ('v3_2', 'file', 'v3.2'), ('v3_2b', 'file', 'v3.2b'),
    ('v3_2b', 'albedo_file', 'v3.2b base reflectance')]):
    p = patch(r, 0.1, version)
    x = 12 + k * 524
    detail.paste(Image.open(source / p[field]).resize((512, 512), Image.Resampling.NEAREST), (x, 65))
    dd.text((x, 580), title, font=font, fill='white')
detail.save(out / 'fiamme_0p1m_magnified.png')
