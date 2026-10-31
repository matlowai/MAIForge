"""Contract checks for the Picea abies calibration fixture.

Proportionate to a prototype: these check that the structure is what it claims to
be and that nothing can quietly stop being a fixture. None of them is evidence of
botanical validity, and passing them says nothing about whether the render looks
like a spruce.
"""
from pathlib import Path
import copy
import math
import sys
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))

import geometry                                          # noqa: E402
import rng                                               # noqa: E402
from evidence import EvidenceError, Ledger, UnavailableMeasurement   # noqa: E402
from fixture import Fixture                              # noqa: E402
from graph import Axis, GraphError                       # noqa: E402

FIXTURE_DIR = ROOT / 'content/visual/benchmarks/picea-abies-calibration-v0'
_CACHE = {}


def built():
    if not _CACHE:
        fixture = Fixture.load(FIXTURE_DIR)
        tree = fixture.build_tree(fixture.solve_envelope()).validate()
        _CACHE.update(fixture=fixture, tree=tree, branch=tree.extract(
            fixture.data['subjects']['branch']['axis']).validate())
    return _CACHE['fixture'], _CACHE['tree'], _CACHE['branch']


class SeedReplayTests(unittest.TestCase):
    def test_rebuilding_from_the_seed_reproduces_the_same_graph(self):
        _, tree, _ = built()
        again = Fixture.load(FIXTURE_DIR)
        rebuilt = again.build_tree(again.solve_envelope())
        self.assertEqual(tree.structural_digest(), rebuilt.structural_digest())
        self.assertEqual([s['id'] for s in tree.sites('ax.w05.p0')],
                         [s['id'] for s in rebuilt.sites('ax.w05.p0')])
        for a, b in zip(tree.sites('ax.w05.p0'), rebuilt.sites('ax.w05.p0')):
            self.assertEqual(a['position'], b['position'])
            self.assertEqual(a['length'], b['length'])

    def test_a_keyed_draw_does_not_depend_on_how_many_draws_came_before(self):
        forwards = [rng.unit(7, f'k{i}') for i in range(6)]
        backwards = [rng.unit(7, f'k{i}') for i in reversed(range(6))]
        self.assertEqual(forwards, list(reversed(backwards)))


class DetailTierTests(unittest.TestCase):
    def test_compiling_every_tier_does_not_mutate_or_regenerate_the_graph(self):
        fixture, _, branch = built()
        before = branch.structural_digest()
        sites_before = branch.sites('ax.w00.p0.s00')
        for tier in geometry.TIERS:
            geometry.compile_tier(branch, tier, fixture.scope, fixture.data, 'branch')
        self.assertEqual(before, branch.structural_digest())
        self.assertIs(sites_before, branch.sites('ax.w00.p0.s00'))

    def test_every_tier_reads_the_same_foliage_sites(self):
        fixture, _, branch = built()
        axis = 'ax.w00.p0.s00'
        near = {s['id']: s['position'] for s in branch.sites(axis)}
        proxies = geometry.compile_tier(branch, 'mid', fixture.scope, fixture.data, 'branch')['proxies']
        for proxy in (p for p in proxies if p['axis'] == axis):
            for site_id in geometry.resolve(branch, proxy):
                self.assertIn(site_id, near)

    def test_mid_and_far_reduce_measured_cost(self):
        fixture, _, branch = built()
        costs = {tier: geometry.compile_tier(branch, tier, fixture.scope, fixture.data, 'branch')['counts']['triangles']
                 for tier in geometry.TIERS}
        self.assertLess(costs['mid'], costs['near'])
        self.assertLess(costs['far'], costs['mid'])

    def test_an_over_budget_tier_refuses_rather_than_silently_reducing(self):
        fixture, tree, _ = built()
        with self.assertRaises(geometry.Budget) as over:
            geometry.compile_tier(tree, 'near', fixture.scope, fixture.data, 'tree')
        self.assertGreater(over.exception.projected, over.exception.limit)
        self.assertEqual(over.exception.detail['sites'], tree.total_site_count())

    def test_an_unknown_tier_fails(self):
        fixture, _, branch = built()
        with self.assertRaises(ValueError):
            geometry.compile_tier(branch, 'ultra', fixture.scope, fixture.data, 'branch')


class ProxyMappingTests(unittest.TestCase):
    def test_every_mid_cluster_resolves_to_sites_it_actually_covers(self):
        fixture, _, branch = built()
        result = geometry.compile_tier(branch, 'mid', fixture.scope, fixture.data, 'branch')
        self.assertTrue(result['proxies'])
        reach = max(fixture.scope.measured('pa.needle.length_m_range'))
        for proxy in result['proxies'][:60]:
            ids = geometry.resolve(branch, proxy)
            self.assertEqual(len(ids), proxy['sites'])
            lookup = {s['id']: s for s in branch.sites(proxy['axis'])}
            for site_id in ids:
                offset = math.dist(lookup[site_id]['position'], proxy['centre'])
                self.assertLessEqual(offset, proxy['half_length'] + proxy['half_width'] + reach)

    def test_every_far_proxy_bounds_the_sites_it_stands_for(self):
        fixture, _, branch = built()
        result = geometry.compile_tier(branch, 'far', fixture.scope, fixture.data, 'branch')
        covered = 0
        for proxy in result['proxies']:
            for site_id in geometry.resolve(branch, proxy):
                axis_id = site_id.split('#')[0]
                site = branch.sites(axis_id)[int(site_id.split('#n')[1])]
                for i in range(3):
                    self.assertLessEqual(abs(site['position'][i] - proxy['centre'][i]),
                                         proxy['half_extents'][i] + 1e-6)
                covered += 1
        self.assertEqual(covered, branch.total_site_count())


class AttachmentTests(unittest.TestCase):
    def test_child_roots_sit_on_their_parent_path_at_the_stored_attachment(self):
        _, tree, _ = built()
        checked = 0
        for axis in tree.axes.values():
            if axis.parent is None:
                continue
            expected, _ = tree.axes[axis.parent].point(axis.parent_arc)
            self.assertLess(math.dist(expected, axis.samples[0][1:4]), 1e-9,
                            f'{axis.id} root drifted off {axis.parent}')
            checked += 1
        self.assertGreater(checked, 1000)

    def test_a_root_moved_off_its_declared_attachment_is_rejected(self):
        _, _, branch = built()
        broken = copy.deepcopy(branch)
        victim = broken.axes['ax.w00.p0.s00']
        drift = list(victim.samples[0])
        drift[2] += 0.05
        victim.samples[0] = tuple(drift)
        with self.assertRaises(GraphError) as failure:
            broken.validate()
        self.assertIn('off its declared attachment', str(failure.exception))


class StructureTests(unittest.TestCase):
    def test_branch_order_and_growth_unit_membership_are_different_things(self):
        _, tree, _ = built()
        for unit in tree.growth_units[:400]:
            axis = tree.axes[unit.axis]
            self.assertIn(axis.order, (0, 1, 2, 3))
            self.assertIsNone(unit.year, 'no chronology was supplied for this population')
        multiple = [a for a in tree.axes.values()
                    if sum(1 for g in tree.growth_units if g.axis == a.id) > 1 and a.order == 1]
        self.assertTrue(multiple, 'a primary should decompose into several growth units')
        self.assertEqual({tree.axes[g.axis].order for g in tree.growth_units if g.axis == multiple[0].id}, {1})

    def test_growth_units_partition_their_axis(self):
        _, _, branch = built()
        for axis_id, axis in branch.axes.items():
            spans = sorted((g.from_s, g.to_s) for g in branch.growth_units if g.axis == axis_id)
            self.assertAlmostEqual(spans[0][0], 0.0)
            self.assertAlmostEqual(spans[-1][1], axis.length, places=6)

    def test_needle_sites_attach_to_twigs_with_finite_dimensions_in_metres(self):
        _, _, branch = built()
        low, high = 0.01, 0.025
        for axis_id in list(branch.axes)[:40]:
            axis = branch.axes[axis_id]
            for site in branch.sites(axis_id):
                self.assertTrue(all(math.isfinite(v) for v in site['position']))
                self.assertGreaterEqual(site['length'], low)
                self.assertLessEqual(site['length'], high)
                self.assertGreater(site['width'], 0)
                centre, radius = axis.point(site['s'])
                self.assertLess(abs(math.dist(centre, site['position']) - radius), 1e-9,
                                'a needle site should sit on its own shoot surface')

    def test_only_twigs_and_distal_primaries_carry_foliage(self):
        _, tree, _ = built()
        self.assertEqual(tree.axes['ax.trunk'].site_count(), 0)
        primary = tree.axes['ax.w00.p0']
        self.assertGreater(primary.foliage['from_s'], 0.0)
        self.assertGreater(tree.axes['ax.w00.p0.s00'].site_count(), 0)


class ExtractionTests(unittest.TestCase):
    def test_the_inspected_branch_is_the_tree_subgraph_not_a_second_generation(self):
        _, tree, branch = built()
        self.assertEqual(branch.extracted_from, 'ax.w00.p0')
        self.assertTrue(set(branch.axes) <= set(tree.axes))
        for axis_id in list(branch.axes)[:25]:
            for local, whole in zip(branch.sites(axis_id), tree.sites(axis_id)):
                self.assertEqual(local['id'], whole['id'])
                for a, b in zip(branch.to_tree_frame(local['position']), whole['position']):
                    self.assertAlmostEqual(a, b, places=9)

    def test_the_extracted_branch_keeps_a_transform_back_to_tree_coordinates(self):
        _, tree, branch = built()
        self.assertIsNotNone(branch.local_to_tree)
        root = tree.axes['ax.w00.p0'].samples[0][1:4]
        self.assertEqual(branch.to_tree_frame((0.0, 0.0, 0.0)), root)


class InvalidGraphTests(unittest.TestCase):
    def _broken(self):
        _, _, branch = built()
        return copy.deepcopy(branch)

    def test_a_missing_parent_fails(self):
        broken = self._broken()
        broken.axes['ax.w00.p0.s00'].parent = 'ax.nowhere'
        self.assertRaises(GraphError, broken.validate)

    def test_a_cycle_fails(self):
        broken = self._broken()
        broken.axes['ax.w00.p0'].parent = 'ax.w00.p0.s00'
        broken.axes['ax.w00.p0'].parent_arc = 0.0
        self.assertRaises(GraphError, broken.validate)

    def test_a_nonfinite_coordinate_fails(self):
        broken = self._broken()
        axis = broken.axes['ax.w00.p0']
        axis.samples[3] = (axis.samples[3][0], float('nan'), 0.0, 0.0, 0.01)
        self.assertRaises(GraphError, broken.validate)

    def test_a_zero_radius_fails_rather_than_passing_as_a_thin_tip(self):
        broken = self._broken()
        axis = broken.axes['ax.w00.p0']
        axis.samples[2] = axis.samples[2][:4] + (0.0,)
        self.assertRaises(GraphError, broken.validate)

    def test_an_unsupported_role_fails(self):
        broken = self._broken()
        broken.axes['ax.w00.p0'].role = 'quaternary'
        self.assertRaises(GraphError, broken.validate)

    def test_a_growth_unit_that_acquires_a_year_fails(self):
        broken = self._broken()
        broken.growth_units[0].year = 2019
        self.assertRaises(GraphError, broken.validate)


class LabelTests(unittest.TestCase):
    def test_a_fixture_cannot_become_a_world_artifact_by_changing_a_label(self):
        _, _, branch = built()
        relabelled = copy.deepcopy(branch)
        relabelled.artifact_class = 'world-plant'
        self.assertRaises(GraphError, relabelled.validate)
        admitted = copy.deepcopy(branch)
        admitted.world_admission = True
        self.assertRaises(GraphError, admitted.validate)

    def test_the_graph_constructor_refuses_a_world_class_outright(self):
        from graph import CalibrationGraph
        with self.assertRaises(GraphError):
            CalibrationGraph('x', 1, {}, None, {}, artifact_class='world-plant')


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.ledger = Ledger.load(FIXTURE_DIR / 'evidence-ledger-v0.json')

    def test_a_declared_gap_raises_instead_of_returning_zero(self):
        for pid in self.ledger.gaps():
            with self.assertRaises(UnavailableMeasurement):
                self.ledger.value(pid)
            self.assertIsNone(self.ledger.row(pid).get('value'))

    def test_a_superseded_choice_cannot_still_be_read(self):
        with self.assertRaises(EvidenceError):
            self.ledger.choice('fx.trunk_base_radius_m')

    def test_a_measurement_cannot_be_read_as_a_choice_or_the_reverse(self):
        with self.assertRaises(EvidenceError):
            self.ledger.choice('pa.tree.height_m')
        with self.assertRaises(EvidenceError):
            self.ledger.measured('fx.whorl_count')

    def test_nothing_claims_to_be_verified_in_a_session_that_could_not_reach_a_source(self):
        self.assertNotIn('verified-measurement', self.ledger.classes())

    def test_reading_a_parameter_the_fixture_never_declared_raises(self):
        fixture, _, _ = built()
        with self.assertRaises(EvidenceError):
            fixture.scope.measured('pa.gu.order_distribution_pct')

    def test_every_generated_structure_carries_provenance(self):
        _, tree, _ = built()
        for axis in tree.axes.values():
            if axis.order:
                self.assertTrue(axis.provenance, f'{axis.id} carries no provenance')
                for pid in axis.provenance:
                    self.assertIn(pid, self.ledger.rows)
        for pid in tree.provenance_of['site']:
            self.assertIn(pid, self.ledger.rows)

    def test_the_reported_angle_range_is_not_used_for_primary_branches(self):
        row = self.ledger.row('pa.branch.nonprimary_angle_deg_range')
        self.assertIn('never a primary insertion angle', row['applicability'])
        self.assertIsNone(self.ledger.rows['pa.branch.primary_insertion_angle_deg'].get('value'))


class EnvelopeTests(unittest.TestCase):
    def test_the_fixture_meets_the_two_envelope_targets_it_declared(self):
        fixture, tree, _ = built()
        tolerance = fixture.scope.choice('fx.envelope_tolerance_m')
        self.assertAlmostEqual(tree.axes['ax.trunk'].length,
                               fixture.scope.measured('pa.tree.height_m'), delta=tolerance)
        self.assertAlmostEqual(fixture.crown_diameter(tree),
                               fixture.scope.measured('pa.tree.crown_diameter_m'), delta=tolerance)

    def test_the_solved_crown_lands_inside_the_reported_population_range(self):
        fixture, tree, _ = built()
        low, high = fixture.ledger.row('pa.tree.crown_diameter_m')['reported_range_m']
        self.assertGreater(fixture.crown_diameter(tree), low)
        self.assertLess(fixture.crown_diameter(tree), high)

    def test_the_trunk_taper_passes_through_the_measured_breast_height_diameter(self):
        fixture, tree, _ = built()
        _, radius = tree.axes['ax.trunk'].point(1.3)
        self.assertAlmostEqual(radius * 2, fixture.scope.measured('pa.tree.dbh_m'), places=3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
