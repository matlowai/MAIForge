import unittest
from copy import deepcopy
from test_growth import config,annual
from model import Growth
from development.growth_delta import growth_delta,invalidation,representation_keys
class DeltaTests(unittest.TestCase):
    def test_real_tick_delta_covers_every_changed_render_key(self):
        g=Growth(config())
        for c in annual()[:9]:g.tick(c)
        before=g.snapshot();g.tick(annual()[9]);after=g.state;d=g.last_delta;inv=invalidation(after,d)
        old,new=representation_keys(before),representation_keys(after)
        changed={k for k in set(old)|set(new) if old.get(k)!=new.get(k)}
        covered={'wood:'+x for x in inv['wood_nodes']}|{'foliage:'+x for x in inv['foliage_nodes']}
        self.assertFalse(changed-covered)
        self.assertTrue(d['added_shoots']);self.assertTrue(d['radius_changed_shoots']);self.assertTrue(d['shed_cohorts'])
        # This young real tick thickens every extant live shape. Do not invent
        # an unchanged-cache claim when biology actually changes all radii.
    def test_single_local_radius_edit_does_not_flush_unrelated_branches(self):
        g=Growth(config())
        for c in annual()[:5]:g.tick(c)
        before=g.snapshot();after=deepcopy(before);after['year']+=1
        sid=next(s for s,u in after['shoots'].items() if after['axes'][u['axis']]['order']==3)
        after['shoots'][sid]['wood_area_m2']*=1.1
        d=growth_delta(before,after);inv=invalidation(after,d)
        self.assertEqual(d['radius_changed_shoots'],[sid]);self.assertLessEqual(len(inv['wood_nodes']),2)
        self.assertEqual(len(inv['far_primary_nodes']),1)
        self.assertLess(len(inv['foliage_nodes']),3)
    def test_delta_refuses_regeneration_and_missing_history(self):
        g=Growth(config());g.tick(annual()[0]);before=g.snapshot();after=deepcopy(before);after['year']+=1
        after['shoots']['trunk/gu001']['end'][2]+=1
        with self.assertRaises(ValueError):growth_delta(before,after)
        after=deepcopy(before);after['year']+=1;after['cohorts'].pop(next(iter(after['cohorts'])))
        with self.assertRaises(ValueError):growth_delta(before,after)
if __name__=='__main__':unittest.main()
