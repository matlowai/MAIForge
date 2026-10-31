"""Build the calibration graph from the fixture and the ledger.

Every botanical number arrives through ScopedLedger, so the fixture file lists
what this builder is allowed to read and the ledger says what kind of thing each
value is. Nothing botanical is written as a literal in here. The literals that do
appear are frame constants, sampling resolution and arithmetic limits.
"""
from pathlib import Path
import json
import math

from evidence import Ledger, ScopedLedger
from graph import Axis, CalibrationGraph
import rng

DEG = math.pi / 180.0
STRATA = ((1, 5), (6, 10), (11, 15), (16, 20))
ANGLE_LIMITS = (12.0, 88.0)     # chosen clamp on any drawn insertion angle, degrees
SAMPLE_STEP_M = 0.05            # centreline sampling resolution, a build setting
MIN_SAMPLES = 3


def _norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) or 1.0
    return (v[0] / n, v[1] / n, v[2] / n)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _axis_frame(tangent):
    """Side and plane-up for an axis: side is horizontal, so left and right lie in
    the roughly horizontal plane a plagiotropic branch develops in."""
    side = _cross(tangent, (0.0, 0.0, 1.0))
    if math.sqrt(sum(c * c for c in side)) < 1e-6:
        side = (1.0, 0.0, 0.0)
    side = _norm(side)
    return side, _norm(_cross(side, tangent))


def _rotate_about(vector, axis, radians):
    c, s = math.cos(radians), math.sin(radians)
    dot = sum(a * b for a, b in zip(vector, axis))
    cr = _cross(axis, vector)
    return tuple(vector[i] * c + cr[i] * s + axis[i] * dot * (1 - c) for i in range(3))


def _stratum_key(from_apex):
    for low, high in STRATA:
        if low <= from_apex <= high:
            return f'{low}-{high}'
    return f'{STRATA[-1][0]}-{STRATA[-1][1]}'


class Fixture:
    def __init__(self, data, ledger, path=None):
        self.data = data
        self.path = path
        self.ledger = ledger
        self.scope = ScopedLedger(ledger, data['consumes'], data['reports']['params'])
        self.notes = {}
        if data.get('artifact_class') != 'visual-calibration-fixture':
            raise ValueError('fixture is not marked a visual calibration fixture')

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        data = json.loads((directory / 'fixture-v0.json').read_text())
        return cls(data, Ledger.load(directory / data['ledger']), directory)

    @property
    def sha256(self):
        from hashlib import sha256
        return sha256(json.dumps(self.data, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

    # -- angle model ------------------------------------------------------
    def _nonprimary_angle(self, key, order):
        """Draw a daughter insertion angle from the measured bin masses.

        The masses are measured; the shape inside a bin is fx.nonprimary_angle_distribution
        and the order shift magnitude is fx.angle_order_shift_deg. Primary branches
        never come through here: the source excluded them from this analysis.
        """
        masses = self.scope.measured('pa.branch.nonprimary_angle_bin_masses')
        low, high = self.scope.measured('pa.branch.nonprimary_angle_deg_range')
        bins = (('below_40', ANGLE_LIMITS[0], low), ('from_40_to_70', low, high), ('at_least_70', high, ANGLE_LIMITS[1]))
        chosen = rng.pick(self.seed, key + '|bin', [b[0] for b in bins], [masses[b[0]] for b in bins])
        span = next(b for b in bins if b[0] == chosen)
        angle = rng.between(self.seed, key + '|within', span[1], span[2])
        angle += self.scope.choice('fx.angle_order_shift_deg') * max(0, order - 2)
        return min(ANGLE_LIMITS[1], max(ANGLE_LIMITS[0], angle))

    def _daughter_direction(self, key, tangent):
        marginals = self.scope.measured('pa.gu.daughter_direction_marginals')
        self.scope.qualitative('pa.crown.branching_is_plagiotropic')
        side, plane_up = _axis_frame(tangent)
        cls = rng.pick(self.seed, key + '|dir', ['left', 'right', 'up'],
                       [marginals['left'], marginals['right'], marginals['up']])
        base = {'left': side, 'right': tuple(-c for c in side), 'up': plane_up}[cls]
        jitter = self.scope.choice('fx.branchlet_azimuth_jitter_deg')
        return cls, _norm(_rotate_about(base, tangent, rng.between(self.seed, key + '|jit', -jitter, jitter) * DEG))

    # -- centreline -------------------------------------------------------
    def _straight_axis(self, origin, direction, length, base_radius):
        taper = self.scope.choice('fx.branch_taper_exponent')
        floor = self.scope.choice('fx.min_wood_radius_m')
        count = max(MIN_SAMPLES, int(math.ceil(length / SAMPLE_STEP_M)) + 1)
        out = []
        for i in range(count):
            s = length * i / (count - 1)
            r = max(floor, base_radius * (1.0 - s / length) ** taper) if s < length else floor
            out.append((s, origin[0] + direction[0] * s, origin[1] + direction[1] * s, origin[2] + direction[2] * s, r))
        return out

    def _drooping_axis(self, origin, azimuth, theta0_deg, length, base_radius):
        """A primary branch that bends progressively downward along its own length."""
        droop = self.scope.choice('fx.primary_droop_deg_per_m')
        taper = self.scope.choice('fx.branch_taper_exponent')
        floor = self.scope.choice('fx.min_wood_radius_m')
        count = max(MIN_SAMPLES, int(math.ceil(length / SAMPLE_STEP_M)) + 1)
        ca, sa = math.cos(azimuth), math.sin(azimuth)
        x, y, z = origin
        out = [(0.0, x, y, z, base_radius)]
        for i in range(1, count):
            s0, s = length * (i - 1) / (count - 1), length * i / (count - 1)
            mid = (s0 + s) * 0.5
            theta = (theta0_deg + droop * mid) * DEG
            step = s - s0
            x += math.sin(theta) * ca * step
            y += math.sin(theta) * sa * step
            z += math.cos(theta) * step
            r = max(floor, base_radius * (1.0 - s / length) ** taper) if s < length else floor
            out.append((s, x, y, z, r))
        return out

    # -- the site rule ----------------------------------------------------
    def site_rule(self, graph, axis, index):
        span = axis.foliage['to_s'] - axis.foliage['from_s']
        count = axis.site_count()
        s = axis.foliage['from_s'] + span * (index + 0.5) / count
        centre, radius = axis.point(s)
        tangent = axis.tangent(s)
        side, plane_up = _axis_frame(tangent)
        azimuth = index * self._needle_step * DEG
        radial = _norm(tuple(side[i] * math.cos(azimuth) + plane_up[i] * math.sin(azimuth) for i in range(3)))
        key = f'{axis.id}|needle|{index}'
        pitch = (self._needle_pitch + rng.between(self.seed, key + '|pitch', -self._pitch_jitter, self._pitch_jitter)) * DEG
        direction = _norm(tuple(tangent[i] * math.cos(pitch) + radial[i] * math.sin(pitch) for i in range(3)))
        low, high = self._needle_range
        return {'id': f'{axis.id}#n{index:04d}', 'axis': axis.id, 'index': index, 's': s,
                'position': tuple(centre[i] + radial[i] * radius for i in range(3)),
                'direction': direction, 'radial': radial,
                'length': rng.between(self.seed, key + '|len', low, high), 'width': self._needle_width}

    # -- the tree ---------------------------------------------------------
    def build_tree(self, crown_scale=None):
        s = self.scope
        self.seed = s.choice('fx.seed')
        height = s.measured('pa.tree.height_m')
        first_whorl = s.choice('fx.first_live_whorl_m')
        adopted = s.measured('pa.tree.lowest_living_whorl_m')
        if abs(first_whorl - adopted) > 1e-9:
            raise ValueError('fx.first_live_whorl_m claims to adopt pa.tree.lowest_living_whorl_m but differs from it')
        leader = s.choice('fx.leader_length_m')
        whorls = s.choice('fx.whorl_count')
        spacing_exp = s.choice('fx.whorl_spacing_exponent')
        taper_exp = s.choice('fx.crown_taper_exponent')
        min_length = s.choice('fx.min_primary_length_m')
        angle_base = s.choice('fx.primary_angle_base_deg')
        angle_apex = s.choice('fx.primary_angle_apex_deg')
        s.qualitative('pa.crown.branch_orientation_by_position')
        s.qualitative('pa.branch.angle_decreases_with_order')
        s.qualitative('pa.gu.axis_decomposes_into_growth_units')
        s.measured('pa.needle.cross_section')
        s.measured('pa.stand.density_stems_per_ha')
        half = s.choice('fx.stratum_radius_from_diameter')
        max_diameter = s.measured('pa.branch.max_diameter_by_stratum_m')
        mean_count = s.measured('pa.branch.mean_count_per_whorl_by_stratum')
        alive = s.measured('pa.branch.alive_fraction_by_stratum')
        within_min = s.choice('fx.within_whorl_radius_min_fraction')
        s.choice('fx.branch_count_realization')
        s.choice('fx.branchlet_direction_scheme')
        s.choice('fx.nonprimary_angle_distribution')
        s.choice('fx.branchlet_azimuth_scheme')
        self._needle_step = s.choice('fx.needle_azimuth_step_deg')
        self._needle_pitch = s.choice('fx.needle_pitch_deg')
        self._pitch_jitter = s.choice('fx.needle_pitch_jitter_deg')
        self._needle_width = s.choice('fx.needle_width_m')
        self._needle_range = s.measured('pa.needle.length_m_range')
        s.choice('fx.needle_length_distribution')
        s.choice('fx.needle_cross_section')
        sites_per_m = s.choice('fx.needle_sites_per_m')
        foliated_from = s.choice('fx.primary_foliated_from_fraction')
        gu_arc = s.choice('fx.growth_unit_arc_m')
        floor = s.choice('fx.min_wood_radius_m')

        graph = CalibrationGraph(self.data['id'], self.seed, self.data['frame'], self.site_rule,
                                 {'site': ('fx.needle_sites_per_m', 'pa.needle.length_m_range', 'fx.needle_width_m',
                                           'fx.needle_pitch_deg', 'fx.needle_azimuth_step_deg', 'pa.needle.cross_section')})
        # trunk: base radius solved so the chosen taper passes through the measured DBH at 1.3 m
        s.choice('fx.trunk_radius_from_dbh')
        dbh, trunk_taper = s.measured('pa.tree.dbh_m'), s.choice('fx.trunk_taper_exponent')
        base_radius = (dbh * half) / ((1.0 - 1.3 / height) ** trunk_taper)
        self.notes['trunk_base_radius_m'] = base_radius
        count = int(math.ceil(height / (SAMPLE_STEP_M * 5))) + 1
        trunk = graph.add_axis(Axis('ax.trunk', 'trunk', None, None,
                                    [(height * i / (count - 1), 0.0, 0.0, height * i / (count - 1),
                                      max(floor, base_radius * (1.0 - i / (count - 1)) ** trunk_taper)) for i in range(count)],
                                    provenance=('pa.tree.height_m', 'pa.tree.dbh_m', 'fx.trunk_radius_from_dbh', 'fx.trunk_taper_exponent')))
        graph.add_growth_units(trunk, gu_arc)

        top = height - leader
        scale = crown_scale if crown_scale is not None else 1.0
        golden = s.choice('fx.whorl_azimuth_step_deg')
        jitter = s.choice('fx.branch_azimuth_jitter_deg')
        secondary_spacing = s.choice('fx.secondary_spacing_m')
        secondary_start = s.choice('fx.secondary_start_m')
        secondary_fraction = s.choice('fx.secondary_length_fraction')
        secondary_max = s.choice('fx.secondary_max_length_m')
        tertiary_spacing = s.choice('fx.tertiary_spacing_m')
        tertiary_start = s.choice('fx.tertiary_start_m')
        tertiary_fraction = s.choice('fx.tertiary_length_fraction')
        tertiary_max = s.choice('fx.tertiary_max_length_m')
        daughter_radius = s.choice('fx.daughter_radius_fraction')

        for w in range(whorls):
            u = w / (whorls - 1)
            z = first_whorl + (top - first_whorl) * (u ** spacing_exp)
            from_apex = whorls - w
            key = _stratum_key(from_apex)
            target = mean_count[key] * alive[key]
            n = int(target) + (1 if rng.bernoulli(self.seed, f'w{w}|count', target - int(target)) else 0)
            length = scale * ((1.0 - u) ** taper_exp) + min_length
            theta = angle_base + (angle_apex - angle_base) * u
            whorl_radius = max_diameter[key] * half
            for b in range(n):
                bkey = f'w{w}|b{b}'
                azimuth = (w * golden + b * 360.0 / n) * DEG + rng.between(self.seed, bkey + '|az', -jitter, jitter) * DEG
                ramp = 1.0 if n == 1 else 1.0 - (1.0 - within_min) * (b / (n - 1))
                primary = graph.add_axis(Axis(
                    f'ax.w{w:02d}.p{b}', 'primary', trunk.id, z,
                    self._drooping_axis((0.0, 0.0, z), azimuth, theta, length, max(floor, whorl_radius * ramp)),
                    foliage={'from_s': length * foliated_from, 'to_s': length, 'sites_per_m': sites_per_m},
                    provenance=('pa.branch.max_diameter_by_stratum_m', 'pa.branch.mean_count_per_whorl_by_stratum',
                                'pa.branch.alive_fraction_by_stratum', 'fx.primary_angle_base_deg',
                                'fx.primary_angle_apex_deg', 'fx.primary_droop_deg_per_m', 'fx.crown_taper_exponent',
                                'pa.crown.branch_orientation_by_position'),
                    stratum=key))
                graph.add_growth_units(primary, gu_arc)
                self._add_daughters(graph, primary, 'secondary', secondary_start, secondary_spacing,
                                    secondary_fraction, secondary_max, daughter_radius, floor, sites_per_m, gu_arc)
        for axis_id in [a for a in list(graph.axes) if graph.axes[a].order == 2]:
            self._add_daughters(graph, graph.axes[axis_id], 'tertiary', tertiary_start, tertiary_spacing,
                                tertiary_fraction, tertiary_max, daughter_radius, floor, sites_per_m, gu_arc)
        return graph

    def _add_daughters(self, graph, parent, role, start, spacing, fraction, cap, radius_fraction, floor, sites_per_m, gu_arc):
        order = {'secondary': 2, 'tertiary': 3}[role]
        index, s = 0, start
        while s < parent.length - spacing * 0.5:
            key = f'{parent.id}|d{index}'
            origin, parent_radius = parent.point(s)
            tangent = parent.tangent(s)
            _, direction = self._daughter_direction(key, tangent)
            angle = self._nonprimary_angle(key, order) * DEG
            heading = _norm(tuple(tangent[i] * math.cos(angle) + direction[i] * math.sin(angle) for i in range(3)))
            length = min(cap, fraction * (parent.length - s))
            if length > SAMPLE_STEP_M:
                axis = graph.add_axis(Axis(
                    f'{parent.id}.{role[0]}{index:02d}', role, parent.id, s,
                    self._straight_axis(origin, heading, length, max(floor, parent_radius * radius_fraction)),
                    foliage={'from_s': 0.0, 'to_s': length, 'sites_per_m': sites_per_m},
                    provenance=('pa.branch.nonprimary_angle_bin_masses', 'pa.gu.daughter_direction_marginals',
                                'pa.crown.branching_is_plagiotropic', 'fx.nonprimary_angle_distribution',
                                'fx.branchlet_direction_scheme', 'fx.angle_order_shift_deg')))
                graph.add_growth_units(axis, gu_arc)
            index, s = index + 1, s + spacing

    # -- envelope solve ---------------------------------------------------
    def crown_diameter(self, graph):
        reach = 0.0
        for axis in graph.axes.values():
            if axis.order == 1:
                for _, x, y, _z, _r in axis.samples:
                    reach = max(reach, math.hypot(x, y))
        return reach * 2.0

    def solve_envelope(self, low=0.5, high=8.0, iterations=34):
        """Solve the primary length scalar so the compiled crown meets its declared
        width target. Fitting to a target is a procedure, not a measurement, and it
        makes no other reported measure agree with its source."""
        target = self.scope.measured('pa.tree.crown_diameter_m')
        for _ in range(iterations):
            mid = (low + high) * 0.5
            if self.crown_diameter(self._skeleton(mid)) < target:
                low = mid
            else:
                high = mid
        return (low + high) * 0.5

    def _skeleton(self, scale):
        """Primary axes only, for the envelope solve: same whorl heights, angles,
        length profile and droop law as the full build, without the daughters and
        foliage that cannot change the crown's outer reach."""
        s = self.scope
        self.seed = s.choice('fx.seed')
        height, first_whorl = s.measured('pa.tree.height_m'), s.choice('fx.first_live_whorl_m')
        leader, whorls = s.choice('fx.leader_length_m'), s.choice('fx.whorl_count')
        spacing_exp, taper_exp = s.choice('fx.whorl_spacing_exponent'), s.choice('fx.crown_taper_exponent')
        min_length = s.choice('fx.min_primary_length_m')
        angle_base, angle_apex = s.choice('fx.primary_angle_base_deg'), s.choice('fx.primary_angle_apex_deg')
        floor, top = s.choice('fx.min_wood_radius_m'), height - leader
        graph = CalibrationGraph(self.data['id'], self.seed, self.data['frame'], self.site_rule, {})
        for w in range(whorls):
            u = w / (whorls - 1)
            z = first_whorl + (top - first_whorl) * (u ** spacing_exp)
            theta = angle_base + (angle_apex - angle_base) * u
            length = scale * ((1.0 - u) ** taper_exp) + min_length
            graph.add_axis(Axis(f'sk.w{w:02d}', 'primary', None, None,
                                self._drooping_axis((0.0, 0.0, z), 0.0, theta, length, max(floor, 0.01))))
        return graph
