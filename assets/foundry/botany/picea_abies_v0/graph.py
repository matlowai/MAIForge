"""The calibration graph: axes, growth units and foliage sites.

Frame is metres, Z up, root collar at the origin. Only what the compilers
actually exercise is represented. There are no buds, carbon budgets or vigour
fields here to make the schema look biological, and there is no calendar: a
growth unit carries year None because no chronology was supplied for this
population.

Branch order and growth unit membership are different things. An axis has one
order and is partitioned into growth units; belonging to a growth unit is not
being a deeper branch.
"""
from hashlib import sha256
import json
import math

ROLE_ORDER = {'trunk': 0, 'primary': 1, 'secondary': 2, 'tertiary': 3}
FIXTURE_CLASS = 'visual-calibration-fixture'
MAX_AXIS_LENGTH_M = 40.0
MAX_RADIUS_M = 1.0


class GraphError(ValueError):
    """The graph is not a thing that can be compiled."""


def _finite(*values):
    return all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)


class Axis:
    __slots__ = ('id', 'role', 'order', 'parent', 'parent_arc', 'samples', 'length',
                 'foliage', 'provenance', 'stratum')

    def __init__(self, id, role, parent, parent_arc, samples, foliage=None, provenance=(), stratum=None):
        self.id = id
        self.role = role
        self.order = ROLE_ORDER[role]
        self.parent = parent
        self.parent_arc = parent_arc
        self.samples = samples          # [(s, x, y, z, radius)]
        self.length = samples[-1][0]
        self.foliage = foliage          # None, or {'from_s', 'to_s', 'sites_per_m'}
        self.provenance = tuple(provenance)
        self.stratum = stratum

    def point(self, s):
        """Centreline position and radius at arc length s, linearly between samples."""
        if s <= 0:
            a = self.samples[0]
            return (a[1], a[2], a[3]), a[4]
        if s >= self.length:
            a = self.samples[-1]
            return (a[1], a[2], a[3]), a[4]
        lo, hi = 0, len(self.samples) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.samples[mid][0] <= s:
                lo = mid
            else:
                hi = mid
        a, b = self.samples[lo], self.samples[hi]
        t = (s - a[0]) / (b[0] - a[0])
        return (a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t, a[3] + (b[3] - a[3]) * t), a[4] + (b[4] - a[4]) * t

    def tangent(self, s):
        lo = min(range(len(self.samples) - 1), key=lambda i: abs(self.samples[i][0] - s))
        lo = min(lo, len(self.samples) - 2)
        a, b = self.samples[lo], self.samples[lo + 1]
        d = (b[1] - a[1], b[2] - a[2], b[3] - a[3])
        n = math.dist((0, 0, 0), d) or 1.0
        return (d[0] / n, d[1] / n, d[2] / n)

    def site_count(self):
        if not self.foliage:
            return 0
        span = self.foliage['to_s'] - self.foliage['from_s']
        return max(0, int(round(span * self.foliage['sites_per_m'])))


class GrowthUnit:
    __slots__ = ('id', 'axis', 'from_s', 'to_s', 'year')

    def __init__(self, id, axis, from_s, to_s):
        self.id = id
        self.axis = axis
        self.from_s = from_s
        self.to_s = to_s
        self.year = None   # unavailable: no chronology was supplied for this population


class CalibrationGraph:
    def __init__(self, fixture_id, seed, frame, site_rule, provenance_of, artifact_class=FIXTURE_CLASS):
        if artifact_class != FIXTURE_CLASS:
            raise GraphError(f'refusing artifact class {artifact_class!r}: this generator only emits {FIXTURE_CLASS}')
        self.artifact_class = artifact_class
        self.world_admission = False
        self.fixture_id = fixture_id
        self.seed = seed
        self.frame = frame
        self.axes = {}
        self.growth_units = []
        self.site_rule = site_rule          # callable(graph, axis, index) -> site dict
        self.provenance_of = provenance_of  # {'site': (...param ids...)}
        self.local_to_tree = None
        self.extracted_from = None
        self._sites = {}

    # -- construction -----------------------------------------------------
    def add_axis(self, axis):
        if axis.id in self.axes:
            raise GraphError(f'duplicate axis id {axis.id}')
        self.axes[axis.id] = axis
        return axis

    def add_growth_units(self, axis, arc):
        cut, index = 0.0, 0
        while cut < axis.length - 1e-9:
            top = min(axis.length, cut + arc)
            self.growth_units.append(GrowthUnit(f'{axis.id}@g{index:02d}', axis.id, cut, top))
            cut, index = top, index + 1
        return index

    # -- foliage ----------------------------------------------------------
    def sites(self, axis_id):
        """Every foliage site on one axis. Pure in (fixture, seed, axis); cached.

        Near, mid and far all reach foliage through this one function, so no tier
        can regrow or perturb what another tier saw.
        """
        cached = self._sites.get(axis_id)
        if cached is None:
            axis = self.axes[axis_id]
            cached = tuple(self.site_rule(self, axis, i) for i in range(axis.site_count()))
            self._sites[axis_id] = cached
        return cached

    def site_ids(self, axis_id, first, last):
        return [s['id'] for s in self.sites(axis_id)[first:last + 1]]

    def total_site_count(self):
        return sum(a.site_count() for a in self.axes.values())

    def foliated_length(self):
        return math.fsum((a.foliage['to_s'] - a.foliage['from_s']) for a in self.axes.values() if a.foliage)

    def axes_by_order(self):
        counts = {}
        for a in self.axes.values():
            counts[a.order] = counts.get(a.order, 0) + 1
        return counts

    def growth_units_by_order(self):
        counts = {}
        for g in self.growth_units:
            order = self.axes[g.axis].order
            counts[order] = counts.get(order, 0) + 1
        return counts

    # -- validation -------------------------------------------------------
    def validate(self, attachment_tolerance=1e-6):
        if self.artifact_class != FIXTURE_CLASS or self.world_admission is not False:
            raise GraphError('a calibration fixture cannot present itself as world state')
        for axis in self.axes.values():
            if axis.role not in ROLE_ORDER:
                raise GraphError(f'{axis.id}: unsupported role {axis.role!r}')
            if len(axis.samples) < 2:
                raise GraphError(f'{axis.id}: an axis needs at least two centreline samples')
            if abs(axis.samples[0][0]) > 1e-12:
                raise GraphError(f'{axis.id}: centreline must start at arc length zero')
            if not (0 < axis.length <= MAX_AXIS_LENGTH_M):
                raise GraphError(f'{axis.id}: invalid axis length {axis.length}')
            previous = -1.0
            for s, x, y, z, r in axis.samples:
                if not _finite(s, x, y, z, r):
                    raise GraphError(f'{axis.id}: nonfinite centreline sample')
                if s <= previous and previous >= 0:
                    raise GraphError(f'{axis.id}: arc length must increase along the centreline')
                if not (0 < r <= MAX_RADIUS_M):
                    raise GraphError(f'{axis.id}: invalid radius {r}; a zero radius is a degenerate tube, not a thin tip')
                previous = s
            if axis.parent is None:
                if axis.order != 0 and axis.id != self.extracted_from:
                    raise GraphError(f'{axis.id}: only the trunk, or the root of an extracted branch, may be parentless')
                continue
            parent = self.axes.get(axis.parent)
            if parent is None:
                raise GraphError(f'{axis.id}: parent {axis.parent} does not exist')
            if axis.order != parent.order + 1:
                raise GraphError(f'{axis.id}: order {axis.order} does not follow parent order {parent.order}')
            if not 0 <= axis.parent_arc <= parent.length + 1e-9:
                raise GraphError(f'{axis.id}: attachment at {axis.parent_arc} lies off parent {parent.id}')
            expected, _ = parent.point(axis.parent_arc)
            root = axis.samples[0][1:4]
            separation = math.dist(expected, root)
            if separation > attachment_tolerance:
                raise GraphError(f'{axis.id}: root sits {separation:.6f} m off its declared attachment '
                                 f'on {parent.id} at arc {axis.parent_arc:.4f}')
        # cycles
        for start in self.axes:
            seen, node, steps = set(), start, 0
            while node is not None:
                if node in seen:
                    raise GraphError(f'cycle in the branch hierarchy at {node}')
                seen.add(node)
                node = self.axes[node].parent
                steps += 1
                if steps > len(self.axes) + 1:
                    raise GraphError('branch hierarchy does not terminate')
        # growth units partition their axis
        spans = {}
        for unit in self.growth_units:
            if unit.axis not in self.axes:
                raise GraphError(f'{unit.id}: growth unit on missing axis {unit.axis}')
            if unit.year is not None:
                raise GraphError(f'{unit.id}: no chronology was supplied for this population')
            spans.setdefault(unit.axis, []).append((unit.from_s, unit.to_s))
        for axis_id, parts in spans.items():
            parts.sort()
            axis = self.axes[axis_id]
            if abs(parts[0][0]) > 1e-9 or abs(parts[-1][1] - axis.length) > 1e-6:
                raise GraphError(f'{axis_id}: growth units do not cover the axis')
            for (_, a), (b, _) in zip(parts, parts[1:]):
                if abs(a - b) > 1e-9:
                    raise GraphError(f'{axis_id}: growth units are not a partition')
        missing = set(self.axes) - set(spans)
        if missing:
            raise GraphError(f'axes with no growth unit decomposition: {sorted(missing)[:3]}')
        return self

    # -- identity ---------------------------------------------------------
    def structural_digest(self):
        h = sha256()
        h.update(json.dumps({'fixture': self.fixture_id, 'seed': self.seed, 'frame': self.frame,
                             'artifact_class': self.artifact_class}, sort_keys=True).encode())
        for axis_id in sorted(self.axes):
            a = self.axes[axis_id]
            h.update(json.dumps([a.id, a.role, a.parent, round(a.parent_arc, 9) if a.parent_arc is not None else None,
                                 [[round(v, 9) for v in s] for s in a.samples],
                                 a.foliage and {k: round(v, 9) if isinstance(v, float) else v for k, v in a.foliage.items()},
                                 list(a.provenance), a.stratum], sort_keys=True).encode())
        for unit in sorted(self.growth_units, key=lambda g: g.id):
            h.update(f'{unit.id}|{unit.axis}|{round(unit.from_s,9)}|{round(unit.to_s,9)}|{unit.year}'.encode())
        return h.hexdigest()

    # -- extraction -------------------------------------------------------
    def extract(self, axis_id):
        """One branch, re-expressed in its own frame, with node ids preserved.

        The inspected branch is this subgraph of the whole tree, not a second
        generation of it: ids, attachment arcs, provenance and foliage sites all
        survive, and local_to_tree carries the geometry back.
        """
        root = self.axes[axis_id]
        origin = root.samples[0][1:4]
        keep, frontier = {axis_id}, [axis_id]
        while frontier:
            current = frontier.pop()
            for other in self.axes.values():
                if other.parent == current and other.id not in keep:
                    keep.add(other.id)
                    frontier.append(other.id)
        out = CalibrationGraph(self.fixture_id, self.seed, dict(self.frame, origin=f'attachment of {axis_id}'),
                               self.site_rule, self.provenance_of)
        out.extracted_from = axis_id
        out.local_to_tree = [[1, 0, 0, origin[0]], [0, 1, 0, origin[1]], [0, 0, 1, origin[2]], [0, 0, 0, 1]]
        shift = lambda s: (s[0], s[1] - origin[0], s[2] - origin[1], s[3] - origin[2], s[4])
        for kept in sorted(keep):
            a = self.axes[kept]
            out.add_axis(Axis(a.id, a.role, a.parent if a.parent in keep else None,
                              a.parent_arc if a.parent in keep else None,
                              [shift(s) for s in a.samples], a.foliage, a.provenance, a.stratum))
        for unit in self.growth_units:
            if unit.axis in keep:
                out.growth_units.append(GrowthUnit(unit.id, unit.axis, unit.from_s, unit.to_s))
        return out

    def to_tree_frame(self, point):
        if self.local_to_tree is None:
            return tuple(point)
        m = self.local_to_tree
        return tuple(m[i][0] * point[0] + m[i][1] * point[1] + m[i][2] * point[2] + m[i][3] for i in range(3))
