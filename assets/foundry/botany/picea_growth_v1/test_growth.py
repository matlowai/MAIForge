import unittest,json
from pathlib import Path
from copy import deepcopy
from model import Growth,digest
from codec import sample_genome
ROOT=Path(__file__).resolve().parents[4]
DATA=ROOT/'content/visual/benchmarks/picea-growth-v1'
def config():return json.loads((DATA/'codec.json').read_text())
def annual():return json.loads((DATA/'annual-inputs.json').read_text())['years']
class AnnualTests(unittest.TestCase):
    def test_replay_prefix_and_birth_records_survive_future_years(self):
        g=Growth(config())
        for c in annual()[:8]:g.tick(c)
        old=g.snapshot();births={k:{n:u[n] for n in ('start','end','length_m','birth_year')} for k,u in old['shoots'].items()}
        for c in annual()[8:12]:g.tick(c)
        for k,u in births.items():self.assertEqual(u,{n:g.state['shoots'][k][n] for n in u})
        self.assertEqual(g.state['tree_id'],old['tree_id']);self.assertEqual(g.state['genome'],old['genome'])
        replay=Growth(config())
        for c in annual()[:8]:replay.tick(c)
        self.assertEqual(digest(replay.state),digest(old))
    def test_cohorts_belong_to_birth_shoot_and_shed_without_deleting_wood(self):
        g=Growth(config())
        for c in annual()[:10]:g.tick(c)
        retired=[co for co in g.state['cohorts'].values() if co['alive_count']==0]
        self.assertTrue(retired)
        for co in retired:
            self.assertIn(co['shoot'],g.state['shoots']);self.assertEqual(co['birth_year'],g.state['shoots'][co['shoot']]['birth_year'])
        self.assertTrue(g.validate())
        self.assertEqual(len(g.state['cohorts']),len(g.state['shoots']))
        self.assertTrue(g.state['cohorts']['trunk/gu010/needles']['alive_count'])
    def test_annual_clusters_follow_variable_leader_segments(self):
        g=Growth(config())
        for c in annual()[:7]:g.tick(c)
        lengths=[r['annual_shoot_length'] for r in g.state['years']];self.assertGreater(max(lengths)-min(lengths),.1)
        for r in g.state['years']:
            for b in r['nodal_branches']:self.assertGreaterEqual(b['fraction'],.8);self.assertLessEqual(b['fraction'],1)
        self.assertTrue(any(r['internodal_branch_count'] for r in g.state['years']))
    def test_low_light_kills_branch_but_keeps_identity_and_wood(self):
        g=Growth(config())
        for c in annual()[:4]:g.tick(c)
        aid=next(a for a in g.state['axes'] if a!='trunk');shoots=g.state['axes'][aid]['shoots'][:]
        for c in annual()[4:8]:g.tick(dict(c,light=.01))
        self.assertFalse(g.state['axes'][aid]['alive'])
        self.assertTrue(all(sid in g.state['shoots'] for sid in shoots))
    def test_bad_conditions_are_transactionally_refused(self):
        g=Growth(config());before=digest(g.state)
        with self.assertRaises(ValueError):g.tick({'year':1})
        self.assertEqual(before,digest(g.state))
        small=config();small['execution_limits']['max_axes']=1;g=Growth(small);before=digest(g.state)
        with self.assertRaises(ValueError):g.tick(annual()[0])
        self.assertEqual(before,digest(g.state))
    def test_incompatible_modules_and_nonfinite_traits_fail(self):
        c=config();c['modules']['foliage']='compound-leaf'
        with self.assertRaises(ValueError):sample_genome(c,1)
        c=config();c['parameters']['branch_sag']['value']=float('nan')
        with self.assertRaises(ValueError):Growth(c)
    def test_outcome_targets_do_not_construct_the_tree(self):
        a=config();b=deepcopy(a);b['outcome_targets'].update(height_m=999,dbh_m=22,crown_diameter_m=123)
        ga,gb=Growth(a),Growth(b)
        for c in annual()[:5]:ga.tick(c);gb.tick(c)
        self.assertEqual(ga.state['axes'],gb.state['axes']);self.assertEqual(ga.state['shoots'],gb.state['shoots'])
    def test_existing_wood_accumulates_and_unknown_genes_fail(self):
        g=Growth(config())
        for c in annual()[:5]:g.tick(c)
        previous={sid:u['wood_area_m2'] for sid,u in g.state['shoots'].items()}
        g.tick(annual()[5])
        self.assertTrue(all(g.state['shoots'][sid]['wood_area_m2']>=area for sid,area in previous.items()))
        bad=config();bad['parameters']['unused_slider']=deepcopy(bad['parameters']['branch_sag'])
        with self.assertRaises(ValueError):Growth(bad)
    def test_execution_capacity_is_not_inherited_genetics(self):
        a=config();b=deepcopy(a);b['execution_limits']['max_axes']=20000
        self.assertEqual(sample_genome(a,7),sample_genome(b,7))
    def test_needle_identity_survives_retirement_and_foliage_has_a_cost(self):
        g=Growth(config());g.tick(annual()[0]);cid=next(iter(g.state['cohorts']));site=g.site(cid,0)
        for c in annual()[1:9]:g.tick(c)
        later=g.site(cid,0);self.assertEqual(site['id'],later['id']);self.assertFalse(later['alive'])
        self.assertGreater(g.state['years'][-1]['resource']['foliage_equivalent_m3'],0)
        with self.assertRaises(ValueError):g.site(cid,999999)
if __name__=='__main__':unittest.main()
