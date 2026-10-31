#!/usr/bin/env python3
"""Export actual generated tile heights and retain failed quality diagnostics."""
import hashlib,json,math,struct,sys
from pathlib import Path
root=Path(sys.argv[1]); block=json.loads((root/'t2/block.json').read_text())
shape=block['block_size_posts']
if isinstance(shape,dict):raise ValueError('Unsupported block shape schema')
nx,ny=map(int,shape);raw=(root/'t2/block_height_mm.i32').read_bytes()
if len(raw)!=4*nx*ny:raise ValueError('Height field length does not match native block shape')
z=struct.unpack('<'+'i'*(nx*ny),raw)
with (root/'terrain.obj').open('w') as f:
 f.write('# MAIForge generated terrain; metres, x east, y north, z up.\n')
 for y in range(ny):
  for x in range(nx):f.write(f'v {2*x} {2*y} {z[y*nx+x]/1000:.3f}\n')
 for y in range(ny-1):
  for x in range(nx-1):
   a=y*nx+x+1;b=a+1;c=a+nx;d=c+1
   f.write(f'f {a} {b} {c}\nf {b} {d} {c}\n')
failed=[]
def walk(v,path):
 if isinstance(v,dict):
  for k,x in v.items():
   if (k=='pass' or k.endswith('_pass')) and x is False:failed.append(path+'/'+k)
   if isinstance(x,(dict,list)):walk(x,path+'/'+k)
 elif isinstance(v,list):
  for i,x in enumerate(v):walk(x,path+'/'+str(i))
for file in ['t1a/checks.json','t1b/w_0_0/checks.json','t2/block.json']:
 walk(json.loads((root/file).read_text()),file)
report={'schema':'maiforge.example-result/v1','dataflow_completed':True,
 'quality_accepted':not failed,'failed_diagnostics':failed,'shape_posts':[nx,ny],
 'height_range_m':[min(z)/1000,max(z)/1000],'height_sha256':hashlib.sha256(raw).hexdigest(),
 'mesh':'terrain.obj','mesh_sha256':hashlib.sha256((root/'terrain.obj').read_bytes()).hexdigest(),
 'scope':'Generated upstream planet/strata -> drainage -> bounded relief -> detail -> geometry export.',
 'limits':['Coarse, short-run interface example; not Earth calibration or settled terrain.',
 'Failed model diagnostics are retained, not reclassified as passing.',
 'No game models, texture packs, game engine or browser viewer is needed.']}
(root/'demo.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if '--strict' in sys.argv and failed:raise SystemExit(2)
