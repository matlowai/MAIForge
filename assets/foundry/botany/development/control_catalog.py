"""Backend-owned control semantics and validation, consumed by schema-driven labs."""
from pathlib import Path
from copy import deepcopy
import json,math
CATALOG=json.loads((Path(__file__).parent/'trait-control-catalog-v1.json').read_text())
MODULES=CATALOG['modules']
PARAMETERS={k for k,r in CATALOG['parameters'].items() if r['semantic_role']=='inherited_trait'}
def check_value(name,value):
    row=CATALOG['parameters'].get(name)
    if row is None:raise ValueError(f'unsupported catalog trait {name}')
    if row['value_type']=='bool':
        if not isinstance(value,bool):raise ValueError(f'{name}: boolean required')
        return value
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError(f'{name}: finite number required')
    if row['value_type']=='int' and value!=int(value):raise ValueError(f'{name}: integer required')
    d=row['domain']
    for bound,cmp in [('min',lambda a,b:a>=b if d.get('min_inclusive',True) else a>b),('max',lambda a,b:a<=b if d.get('max_inclusive',True) else a<b)]:
        if d.get(bound) is not None and not cmp(value,d[bound]):raise ValueError(f'{name}: outside catalog {bound}')
    return value

def check_constraints(genes):
    for row in CATALOG['constraints']:
        value=sum(genes[k]*v for k,v in row['terms'].items())
        if 'min' in row and (value<row['min'] or value==row['min'] and not row.get('min_inclusive',True)):raise ValueError(row['id'])
        if 'max' in row and (value>row['max'] or value==row['max'] and not row.get('max_inclusive',True)):raise ValueError(row['id'])

def catalog_for(codec):
    if codec['modules']!=MODULES:raise ValueError('unsupported developmental modules')
    result=deepcopy(CATALOG)
    for name in PARAMETERS:
        source=codec['parameters'][name]
        result['parameters'][name].update(evidence_class=source['evidence_class'],provenance=source['provenance'])
    return result
