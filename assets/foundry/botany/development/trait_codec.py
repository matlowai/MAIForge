"""Minimal trait-program boundary, deliberately not a universal plant grammar."""
from hashlib import sha256
import json,math
from .control_catalog import CATALOG,MODULES,PARAMETERS,check_value,check_constraints
def sample_genome(codec,seed):
    if set(codec['parameters'])!=PARAMETERS:raise ValueError('missing or unsupported trait parameter')
    if codec['modules']!=MODULES:raise ValueError('unsupported or incompatible developmental modules')
    genes={}
    for name,row in codec['parameters'].items():
        if 'distribution' in row:
            d=row['distribution']
            if d['kind'] not in CATALOG['parameters'][name]['supported_distribution_kinds']:raise ValueError('unsupported inherited distribution')
            check_value(name,d['low']);check_value(name,d['high'])
            if d['low']>d['high']:raise ValueError('reversed inherited distribution')
            u=int.from_bytes(sha256(f'genome|{seed}|{name}'.encode()).digest()[:8],'big')/2**64
            value=d['low']+(d['high']-d['low'])*u
        else:value=row['value']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('invalid inherited trait')
        check_value(name,value)
        genes[name]=value
    check_constraints(genes)
    cap=codec['execution_limits']['max_axes']
    if not isinstance(cap,int) or isinstance(cap,bool):raise ValueError('integer execution capacity required')
    check_value('execution.max_axes',cap)
    return {'codec_id':codec['id'],'taxon':codec['taxon'],'inheritance_seed':seed,'traits':genes,'modules':dict(MODULES),
            'id':sha256(json.dumps(genes,sort_keys=True).encode()).hexdigest(),
            'evidence':'experimental inherited choices; not a fitted species distribution'}
