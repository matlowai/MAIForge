"""Derived growth receipts and representation dependency invalidation.

Receipts are not part of organism identity. They describe two immutable states;
no renderer tier, cache state, or camera feeds back into development.
"""
from hashlib import sha256
import json

def fingerprint(value):return sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def growth_delta(before,after):
    if before['tree_id']!=after['tree_id'] or before['genome']['id']!=after['genome']['id']:
        raise ValueError('a growth delta cannot change organism/genome identity')
    if after['year']!=before['year']+1:raise ValueError('growth delta requires adjacent annual states')
    if before['world_admission'] is not False or after['world_admission'] is not False:raise ValueError('calibration delta cannot consume world state')
    for collection in ('axes','shoots','cohorts'):
        if set(before[collection])-set(after[collection]):raise ValueError('historical objects cannot disappear')
    added=lambda name:sorted(set(after[name])-set(before[name]))
    common=lambda name:sorted(set(before[name])&set(after[name]))
    radius=[]
    for sid in common('shoots'):
        a,b=before['shoots'][sid],after['shoots'][sid]
        if any(a[k]!=b[k] for k in ('id','axis','birth_year','start','end','length_m','birth_area_m2')):
            raise ValueError('existing shoot birth geometry was morphed')
        if b['wood_area_m2']<a['wood_area_m2']:raise ValueError('accumulated wood area shrank')
        if a['wood_area_m2']!=b['wood_area_m2']:radius.append(sid)
    population=[]
    for cid in common('cohorts'):
        a,b=before['cohorts'][cid],after['cohorts'][cid]
        if any(a[k]!=b[k] for k in ('id','shoot','axis','birth_year','birth_count')):raise ValueError('cohort birth identity changed')
        if a['alive_count']!=b['alive_count']:population.append(cid)
    return {'schema':'growth-delta/v1','tree_id':after['tree_id'],'genome_id':after['genome']['id'],
            'from_year':before['year'],'to_year':after['year'],'world_admission':False,
            'added_axes':added('axes'),'changed_axes':[aid for aid in common('axes') if before['axes'][aid]!=after['axes'][aid]],
            'dead_axes':[aid for aid in common('axes') if before['axes'][aid]['alive'] and not after['axes'][aid]['alive']],
            'added_shoots':added('shoots'),'radius_changed_shoots':radius,'added_cohorts':added('cohorts'),
            'population_changed_cohorts':population,'shed_cohorts':[cid for cid in population if after['cohorts'][cid]['alive_count']==0]}

def invalidation(state,delta):
    if state['tree_id']!=delta['tree_id'] or state['year']!=delta['to_year'] or state['genome']['id']!=delta['genome_id']:raise ValueError('invalidation revision mismatch')
    changed_radius=set(delta['radius_changed_shoots']);wood=set(delta['added_shoots'])|changed_radius
    for aid in delta['dead_axes']:wood.update(state['axes'][aid]['shoots'])
    # A woody segment's terminal radius comes from the NEXT annual segment.
    # Appending/thickening that successor invalidates the predecessor's taper too.
    predecessors={b:a for axis in state['axes'].values() for a,b in zip(axis['shoots'],axis['shoots'][1:])}
    wood.update(predecessors[sid] for sid in list(wood) if sid in predecessors)
    foliage=set(delta['added_cohorts'])|set(delta['population_changed_cohorts'])
    foliage.update(cid for cid,c in state['cohorts'].items() if c['shoot'] in changed_radius and c['alive_count'])
    bounds=set();primary=set()
    for sid in wood:
        aid=state['shoots'][sid]['axis']
        while aid is not None:
            bounds.add(aid);axis=state['axes'][aid]
            if axis['order']==1:primary.add(aid)
            aid=axis['parent']
    for cid in foliage:
        aid=state['cohorts'][cid]['axis']
        while aid is not None:
            bounds.add(aid);axis=state['axes'][aid]
            if axis['order']==1:primary.add(aid)
            aid=axis['parent']
    return {'schema':'representation-invalidation/v1','wood_nodes':sorted(wood),'foliage_nodes':sorted(foliage),
            'far_primary_nodes':sorted(primary),'bounds_changed_axes':sorted(bounds)}

def representation_keys(state):
    keys={};successor={a:b for axis in state['axes'].values() for a,b in zip(axis['shoots'],axis['shoots'][1:])}
    for sid,u in state['shoots'].items():
        end_area=state['shoots'][successor[sid]]['wood_area_m2'] if sid in successor else u['birth_area_m2']
        keys['wood:'+sid]=fingerprint(['wood/v1',u['start'],u['end'],u['wood_area_m2'],end_area,state['axes'][u['axis']]['alive']])
    traits=state['genome']['traits'];appearance=[traits[k] for k in ('needle_length_m','needle_width_m','needle_pitch_deg','needle_divergence_deg')]
    for cid,c in state['cohorts'].items():
        if not c['alive_count']:continue
        u=state['shoots'][c['shoot']]
        keys['foliage:'+cid]=fingerprint(['foliage/v1',c['birth_count'],c['alive_count'],u['start'],u['end'],u['wood_area_m2'],appearance])
    return keys
