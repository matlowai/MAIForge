"""Persistent annual-development hypothesis; never world state or fitted botany.

Birth geometry, shoot identities and cohort populations are append-only. Annual
wood area and live/dead status are explicit state revisions with event receipts.
The finished target height/DBH are not inputs to this model.
"""
from copy import deepcopy
from hashlib import sha256
import json, math
from .trait_codec import sample_genome
from .growth_delta import growth_delta


def draw(seed,key):
    return int.from_bytes(sha256(f'{seed}|{key}'.encode()).digest()[:8],'big')/2**64

def unit(v):
    n=math.sqrt(sum(x*x for x in v))
    if n<=0:raise ValueError('zero direction')
    return [x/n for x in v]

def digest(value):return sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


class WoodyAnnualBackend:
    def __init__(self,config):
        self.config=deepcopy(config)
        if config['artifact_class']!='botanical-growth-calibration' or config['world_admission'] is not False:
            raise ValueError('growth calibration is not world input')
        self.genome=sample_genome(config,config['seed']);self.p=self.genome['traits']
        for k,r in config['parameters'].items():
            if not r.get('provenance') or r['evidence_class'] not in ('explicit-calibration-choice','supplied-unverified-claim'):
                raise ValueError(f'{k}: missing or unsupported evidence')
            if 'value' in r and (not isinstance(r['value'],(int,float)) or not math.isfinite(r['value'])):raise ValueError(f'{k}: nonfinite choice')
        self.state={'schema':'woody-annual-state/v1','artifact_class':config['artifact_class'],'world_admission':False,
                    'tree_id':config['tree_id'],'genome':deepcopy(self.genome),'seed':config['seed'],'year':0,'axes':{},'shoots':{},'cohorts':{},
                    'events':[],'years':[],'root_allocation_m3':0.0,'config_digest':digest(config)}
        self._axis('trunk',None,None,0,0,[0,0,1],'leader')
        self.last_delta=None

    def _axis(self,aid,parent,attachment,order,year,direction,kind):
        if aid in self.state['axes']:raise ValueError('duplicate axis')
        a={'id':aid,'parent':parent,'attachment':attachment,'order':order,'birth_year':year,'direction':unit(direction),
           'kind':kind,'alive':True,'death_year':None,'low_light_years':0,'shoots':[],'light':1.,'vigor':1.,
           'provenance':['branch_angle_deg','terminal_upturn','branch_sag','bud_probability']}
        self.state['axes'][aid]=a
        return a

    def _tip(self,axis):
        if axis['shoots']:return self.state['shoots'][axis['shoots'][-1]]['end'][:]
        if axis['parent'] is None:return [0.,0.,0.]
        at=axis['attachment'];shoot=self.state['shoots'][at['shoot']]
        return [a+(b-a)*at['fraction'] for a,b in zip(shoot['start'],shoot['end'])]

    def _shoot(self,axis,year,length,birth_radius):
        sid=f"{axis['id']}/gu{year:03d}"
        if sid in self.state['shoots']:raise ValueError('shoot already born this year')
        start=self._tip(axis);direction=axis['direction'][:]
        if axis['order']:
            age=year-axis['birth_year'];sag=self.p['branch_sag']*min(age,self.p['sag_maturity_years'])
            # New terminal growth responds to age and a terminal bias; old birth
            # geometry is not morphed. Mechanical aging of existing wood is future.
            direction=unit([direction[0],direction[1],direction[2]-sag+self.p['terminal_upturn']])
        end=[start[i]+direction[i]*length for i in range(3)]
        area=math.pi*birth_radius**2
        shoot={'id':sid,'axis':axis['id'],'birth_year':year,'start':start,'end':end,'length_m':length,
               'birth_area_m2':area,'wood_area_m2':area,'required_pipe_m2':0.,'active_pipe_m2':0.,'disused_pipe_m2':0.,'rings':[],
               'provenance':['leader_increment_m','branch_increment_m','shoot_birth_radius_m']}
        self.state['shoots'][sid]=shoot;axis['shoots'].append(sid)
        self.state['events'].append({'year':year,'type':'shoot_added','id':sid,'birth_digest':digest({k:shoot[k] for k in ('id','axis','birth_year','start','end','length_m','birth_area_m2')})})
        cid=sid+'/needles';count=max(1,round(length*self.p['needle_density_per_m']))
        cohort={'id':cid,'shoot':sid,'axis':axis['id'],'birth_year':year,'birth_count':count,'alive_count':count,'shed_year':None,
                'provenance':['needle_density_per_m','needle_retention_years','needle_length_m','needle_width_m','needle_divergence_deg','needle_pitch_deg']}
        self.state['cohorts'][cid]=cohort
        return shoot

    def tick(self,conditions):
        p=self.p;year=self.state['year']+1
        required={'year','light','extension_factor','new_wood_m3','competition'}
        if set(conditions)!=required or conditions['year']!=year:raise ValueError('missing, extra or out-of-sequence annual conditions')
        if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in conditions.values()):raise ValueError('invalid annual conditions')
        if not 0<=conditions['light']<=1 or not 0<=conditions['competition']<=1 or conditions['extension_factor']<0 or conditions['new_wood_m3']<=0:raise ValueError('invalid annual resource/light input')
        # Work transactionally. A refused allocation leaves the previous year intact.
        before=deepcopy(self.state)
        try:
            self._tick(year,conditions)
            self.last_delta=growth_delta(before,self.state)
        except Exception:self.state=before;raise
        return self.state

    def _tick(self,year,c):
        p=self.p;s=self.state;seed=s['seed'];height=self._tip(s['axes']['trunk'])[2]
        for axis in sorted(s['axes'].values(),key=lambda a:a['order']):
            z=self._tip(axis)[2];depth=max(0,height-z)
            light=c['light']*math.exp(-p['shade_per_m']*depth*(p['shade_baseline']+c['competition']))
            axis['light']=light;axis['vigor']=light*c['extension_factor']
            if not axis['alive']:continue
            axis['low_light_years']=axis['low_light_years']+1 if light<p['death_light_threshold'] else 0
            parent_dead=axis['parent'] is not None and not s['axes'][axis['parent']]['alive']
            if axis['order'] and (parent_dead or axis['low_light_years']>=p['death_delay_years']):
                axis['alive']=False;axis['death_year']=year;s['events'].append({'year':year,'type':'axis_died','id':axis['id'],'cause':'parent_death' if parent_dead else 'low_light_fixture_rule'})
        for cohort in s['cohorts'].values():
            if cohort['alive_count'] and (year-cohort['birth_year']>=p['needle_retention_years'] or not s['axes'][cohort['axis']]['alive']):
                cohort['alive_count']=0;cohort['shed_year']=year;s['events'].append({'year':year,'type':'cohort_shed','id':cohort['id']})
        # Allocation shares are display-model resource choices in wood-volume units;
        # this is NOT a carbon-assimilation equation or a calibrated physiological model.
        available=c['new_wood_m3'];roots=available*p['root_allocation_fraction'];extension_budget=available*p['extension_allocation_fraction']
        birth_radius=p['shoot_birth_radius_m'];area=math.pi*birth_radius**2
        existing=list(s['axes'].values());proposals=[]
        for axis in existing:
            if not axis['alive']:continue
            age=year-axis['birth_year']
            if axis['order']==0:
                length=p['leader_increment_m']*(1-math.exp(-year/p['juvenile_ramp_years']))*c['extension_factor']
            else:
                length=p['branch_increment_m']*p['order_extension_factor']**(axis['order']-1)*math.exp(-age/p['branch_extension_decay_years'])*axis['vigor']
            length*=p['extension_noise_low']+(p['extension_noise_high']-p['extension_noise_low'])*draw(seed,f"{axis['id']}|extension|{year}")
            if length>p['minimum_shoot_m']:proposals.append((axis,length))
        demand=sum(length*area+round(length*p['needle_density_per_m'])*p['foliage_cost_equivalent_m3_per_needle'] for a,length in proposals)
        scale=min(1,extension_budget/demand) if demand else 1
        new_shoots=[];spent=0
        for axis,length in proposals:
            if length*scale<=p['minimum_shoot_m']:continue
            shoot=self._shoot(axis,year,length*scale,birth_radius);new_shoots.append(shoot);spent+=shoot['length_m']*area
        leader=next((u for u in new_shoots if u['axis']=='trunk'),None)
        nodal=[];internodal=[];born_axes=[]
        if leader:
            n=int(p['nodal_branch_mean'])+(draw(seed,f'nodal-count|{year}')<p['nodal_branch_mean']%1)
            phase=year*p['year_phase_deg']*math.pi/180
            for j in range(n):
                fraction=1-p['nodal_span_fraction']*draw(seed,f'nodal-height|{year}|{j}')
                az=phase+math.tau*(j+(2*p['azimuth_jitter_fraction']*draw(seed,f'az|{year}|{j}')-p['azimuth_jitter_fraction']))/n
                angle=p['branch_angle_deg']*math.pi/180
                aid=f"primary.y{year:03d}.n{j:02d}"
                a=self._axis(aid,'trunk',{'shoot':leader['id'],'fraction':fraction},1,year,[math.cos(az)*math.sin(angle),math.sin(az)*math.sin(angle),math.cos(angle)],'nodal')
                born_axes.append(a);nodal.append({'axis':aid,'z':leader['start'][2]+leader['length_m']*fraction,'fraction':fraction})
            n_internodal=sum(draw(seed,f'internode|{year}|{j}')<p['internodal_probability'] for j in range(int(p['internodal_trials'])))
            for j in range(n_internodal):
                fraction=p['internodal_min_fraction']+(1-p['nodal_span_fraction']-2*p['internodal_min_fraction'])*draw(seed,f'ih|{year}|{j}');az=math.tau*draw(seed,f'ia|{year}|{j}')
                aid=f'primary.y{year:03d}.i{j:02d}'
                a=self._axis(aid,'trunk',{'shoot':leader['id'],'fraction':fraction},1,year,[math.cos(az),math.sin(az),p['internodal_vertical_bias']],'internodal')
                born_axes.append(a);internodal.append({'axis':aid,'fraction':fraction,'z':leader['start'][2]+leader['length_m']*fraction})
        for shoot in new_shoots:
            parent=s['axes'][shoot['axis']]
            if not 0<parent['order']<p['maximum_axis_order']:continue
            for j in range(int(p['daughter_trials'])):
                if draw(seed,f"{shoot['id']}|bud|{j}")>=p['bud_probability']*parent['light']:continue
                aid=f"{parent['id']}.y{year:03d}.d{j}"
                d=parent['direction'];az=math.atan2(d[1],d[0])+(1 if j%2 else -1)*p['daughter_angle_deg']*math.pi/180
                born_axes.append(self._axis(aid,parent['id'],{'shoot':shoot['id'],'fraction':p['daughter_attachment_min']+(p['daughter_attachment_max']-p['daughter_attachment_min'])*draw(seed,aid+'|at')},parent['order']+1,year,[math.cos(az),math.sin(az),p['daughter_vertical_bias']],'daughter'))
        if len(s['axes'])>self.config['execution_limits']['max_axes']:raise ValueError('development exceeds explicit axis capacity; no biological pruning applied')
        for axis in born_axes:
            length=p['new_branch_shoot_m']*p['order_extension_factor']**(axis['order']-1)*c['extension_factor']
            cost=length*area
            if spent+cost>extension_budget:raise ValueError('new bud development exceeds extension allocation')
            self._shoot(axis,year,length,birth_radius);spent+=cost
        # Supported foliage flows through actual attachment shoots and the stem
        # above each segment. Dead wood persists even when active pipe demand falls.
        foliage_cost=sum(co['birth_count']*p['foliage_cost_equivalent_m3_per_needle'] for co in s['cohorts'].values() if co['birth_year']==year)
        if spent+foliage_cost>extension_budget:raise ValueError('new foliage exceeds development allocation')
        own={sid:0 for sid in s['shoots']};attached={sid:[] for sid in s['shoots']}
        for co in s['cohorts'].values():own[co['shoot']]+=co['alive_count']
        for a in s['axes'].values():
            if a['parent'] is not None:attached[a['attachment']['shoot']].append(a['id'])
        support={};axis_support={}
        for axis in sorted(s['axes'].values(),key=lambda a:-a['order']):
            above=0
            for sid in reversed(axis['shoots']):
                above+=own[sid]+sum(axis_support[aid] for aid in attached[sid]);support[sid]=above
            axis_support[axis['id']]=above
        deficits={};radial_budget=available-roots-spent-foliage_cost
        for sid,shoot in s['shoots'].items():
            active=support[sid]*p['pipe_area_per_needle_m2'];shoot['required_pipe_m2']=active
            deficits[sid]=max(0,active-shoot['wood_area_m2'])*shoot['length_m']
        total=sum(deficits.values());allocation=min(1,radial_budget/total) if total else 0;radial=0
        for sid,shoot in s['shoots'].items():
            volume=deficits[sid]*allocation;inc=volume/shoot['length_m'];shoot['wood_area_m2']+=inc
            shoot['active_pipe_m2']=min(shoot['wood_area_m2'],shoot['required_pipe_m2'])
            shoot['disused_pipe_m2']=shoot['wood_area_m2']-shoot['active_pipe_m2']
            if inc:shoot['rings'].append({'year':year,'area_increment_m2':inc,'supported_needles':support[sid]})
            radial+=volume
        s['root_allocation_m3']+=roots;s['year']=year
        report={'year':year,'conditions':deepcopy(c),'annual_shoot_length':leader['length_m'] if leader else 0,
                'nodal_cluster_span':leader['length_m']*p['nodal_span_fraction'] if leader else 0,
                'nodal_branch_count':len(nodal),'internodal_branch_count':len(internodal),'nodal_branches':nodal,'internodal_branches':internodal,
                'axes_added':[a['id'] for a in born_axes],'resource':{'available_m3':available,'root_m3':roots,'extension_m3':spent,'foliage_equivalent_m3':foliage_cost,'radial_m3':radial,'unallocated_m3':max(0,available-roots-spent-foliage_cost-radial)},
                'pipe_shortfall_m3':max(0,total-radial),'geometry_aging':'new terminal segments only; existing load-induced bending not simulated'}
        s['years'].append(report);s['events'].append({'year':year,'type':'annual_tick','conditions_digest':digest(c),'new_axes':len(born_axes)})
        self.validate()

    def validate(self):
        s=self.state
        for axis in s['axes'].values():
            if axis['parent'] is not None:
                parent=s['axes'][axis['parent']];at=axis['attachment'];u=s['shoots'][at['shoot']]
                assert u['axis']==parent['id'] and axis['birth_year']>=u['birth_year']
                assert 0<=at['fraction']<=1
            previous=None
            for sid in axis['shoots']:
                shoot=s['shoots'][sid]
                if previous:assert shoot['start']==previous['end'] and shoot['birth_year']>previous['birth_year']
                assert shoot['wood_area_m2']>=shoot['birth_area_m2']>0
                assert abs(shoot['active_pipe_m2']+shoot['disused_pipe_m2']-shoot['wood_area_m2'])<1e-12
                previous=shoot
        for co in s['cohorts'].values():
            assert co['birth_year']==s['shoots'][co['shoot']]['birth_year']
            assert co['alive_count'] in (0,co['birth_count'])
        for r in s['years']:
            budget=r['resource'];assert abs(budget['available_m3']-sum(budget[k] for k in ('root_m3','extension_m3','foliage_equivalent_m3','radial_m3','unallocated_m3')))<1e-10
        return True

    def site(self,cohort_id,index):
        co=self.state['cohorts'][cohort_id]
        if not isinstance(index,int) or not 0<=index<co['birth_count']:raise ValueError('needle index outside cohort')
        u=self.state['shoots'][co['shoot']];t=(index+.5)/co['birth_count']
        d=unit([b-a for a,b in zip(u['start'],u['end'])]);ref=[0,0,1] if abs(d[2])<.9 else [1,0,0]
        cross=lambda a,b:[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
        side=unit(cross(d,ref));up=unit(cross(side,d));angle=math.radians(index*self.p['needle_divergence_deg']);pitch=math.radians(self.p['needle_pitch_deg'])
        radial=[side[i]*math.cos(angle)+up[i]*math.sin(angle) for i in range(3)];radius=math.sqrt(u['wood_area_m2']/math.pi)
        position=[u['start'][i]+(u['end'][i]-u['start'][i])*t+radial[i]*radius for i in range(3)]
        return {'id':f'{cohort_id}#n{index:05d}','cohort':cohort_id,'axis':co['axis'],'birth_year':co['birth_year'],
                'alive':bool(co['alive_count']),'position':position,'direction':unit([d[i]*math.cos(pitch)+radial[i]*math.sin(pitch) for i in range(3)]),
                'length':self.p['needle_length_m'],'width':self.p['needle_width_m']}

    def snapshot(self):return deepcopy(self.state)
