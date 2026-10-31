"""Near, mid and far compilers over one calibration graph.

All three read the SAME graph and the same foliage sites. None of them regrows,
reseeds or perturbs anything: mid and far reduce, and every reduction keeps a
mapping back to the source sites it stands for.

Near builds real woody tubes and individual needles. Mid keeps the woody topology
and replaces each short run of needles with one bounded twig scale cluster. Far
keeps the trunk and primaries and aggregates the crown's occupied regions into
box proxies. Budgets are explicit test settings; a tier that exceeds one is
reported unavailable with its measured projection rather than quietly reduced.
"""
import math

TIERS = ('near', 'mid', 'far')


class Budget(Exception):
    def __init__(self, tier, projected, limit, detail):
        super().__init__(f'{tier} needs {projected} triangles, budget is {limit}')
        self.tier, self.projected, self.limit, self.detail = tier, projected, limit, detail


class Mesh:
    __slots__ = ('name', 'role', 'positions', 'normals', 'indices')

    def __init__(self, name, role):
        self.name, self.role = name, role
        self.positions, self.normals, self.indices = [], [], []

    def vertex(self, p, n):
        i = len(self.positions) // 3
        self.positions.extend(p)
        self.normals.extend(n)
        return i

    def triangle(self, a, b, c):
        self.indices.extend((a, b, c))

    @property
    def vertex_count(self):
        return len(self.positions) // 3

    @property
    def triangle_count(self):
        return len(self.indices) // 3

    def payload(self, digits=5):
        return {'name': self.name, 'role': self.role,
                'positions': [round(v, digits) for v in self.positions],
                'normals': [round(v, 5) for v in self.normals],
                'indices': self.indices,
                'vertices': self.vertex_count, 'triangles': self.triangle_count}


def _norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) or 1.0
    return (v[0] / n, v[1] / n, v[2] / n)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _perpendicular(t):
    other = (0.0, 0.0, 1.0) if abs(t[2]) < 0.9 else (1.0, 0.0, 0.0)
    return _norm(_cross(t, other))


def _frames(samples):
    """Rotation minimising frame along a centreline, so a tube does not twist."""
    points = [(s[1], s[2], s[3]) for s in samples]
    tangents = []
    for i in range(len(points)):
        a = points[max(0, i - 1)]
        b = points[min(len(points) - 1, i + 1)]
        tangents.append(_norm((b[0] - a[0], b[1] - a[1], b[2] - a[2])))
    normals = [_perpendicular(tangents[0])]
    for i in range(1, len(points)):
        previous, t = normals[-1], tangents[i]
        projected = tuple(previous[j] - t[j] * sum(previous[k] * t[k] for k in range(3)) for j in range(3))
        if math.sqrt(sum(c * c for c in projected)) < 1e-9:
            projected = _perpendicular(t)
        normals.append(_norm(projected))
    return points, tangents, normals


def tube(mesh, axis, segments):
    points, tangents, normals = _frames(axis.samples)
    rings = []
    for i, (point, tangent, normal) in enumerate(zip(points, tangents, normals)):
        binormal = _cross(tangent, normal)
        radius = axis.samples[i][4]
        ring = []
        for k in range(segments):
            a = 2 * math.pi * k / segments
            direction = _norm(tuple(normal[j] * math.cos(a) + binormal[j] * math.sin(a) for j in range(3)))
            ring.append(mesh.vertex(tuple(point[j] + direction[j] * radius for j in range(3)), direction))
        rings.append(ring)
    for lower, upper in zip(rings, rings[1:]):
        for k in range(segments):
            n = (k + 1) % segments
            mesh.triangle(lower[k], lower[n], upper[k])
            mesh.triangle(lower[n], upper[n], upper[k])
    return len(rings)


def needle(mesh, site, sides=4):
    """One needle: a tapered prism from a base ring on the wood to a single tip.

    Four sided because the flora describes a four angled needle; the proportions
    are a calibration choice. There is no base cap, since the base sits on wood.
    """
    direction = site['direction']
    base = site['position']
    half = site['width'] * 0.5
    u = _perpendicular(direction)
    v = _cross(direction, u)
    tip = mesh.vertex(tuple(base[j] + direction[j] * site['length'] for j in range(3)), direction)
    ring = []
    for k in range(sides):
        a = 2 * math.pi * k / sides
        out = _norm(tuple(u[j] * math.cos(a) + v[j] * math.sin(a) for j in range(3)))
        ring.append(mesh.vertex(tuple(base[j] + out[j] * half for j in range(3)), out))
    for k in range(sides):
        mesh.triangle(ring[k], ring[(k + 1) % sides], tip)
    return sides + 1


def buckets(graph, arc, needle_low, needle_high, pitch_deg):
    """Twig scale clusters: a short run of consecutive sites on one axis.

    A bucket is addressed by axis and site index range, so the sites it stands
    for can always be resolved back out of the graph.
    """
    mean = (needle_low + needle_high) * 0.5
    lateral = mean * math.sin(pitch_deg * math.pi / 180.0)
    forward = mean * math.cos(pitch_deg * math.pi / 180.0)
    out = []
    for axis_id in sorted(graph.axes):
        axis = graph.axes[axis_id]
        count = axis.site_count()
        if not count:
            continue
        span = axis.foliage['to_s'] - axis.foliage['from_s']
        groups = max(1, int(math.ceil(span / arc)))
        per = count / groups
        for g in range(groups):
            first, last = int(g * per), min(count - 1, int((g + 1) * per) - 1)
            if last < first:
                continue
            s0 = axis.foliage['from_s'] + span * g / groups
            s1 = axis.foliage['from_s'] + span * (g + 1) / groups
            centre, radius = axis.point((s0 + s1) * 0.5)
            out.append({'id': f'{axis_id}~b{g:03d}', 'axis': axis_id, 'first': first, 'last': last,
                        'sites': last - first + 1, 'from_s': s0, 'to_s': s1,
                        'centre': centre, 'tangent': axis.tangent((s0 + s1) * 0.5),
                        'half_length': (s1 - s0) * 0.5 + forward, 'half_width': radius + lateral})
    return out


def cards(mesh, bucket, count):
    t = bucket['tangent']
    u = _perpendicular(t)
    v = _cross(t, u)
    for k in range(count):
        a = math.pi * k / count
        side = _norm(tuple(u[j] * math.cos(a) + v[j] * math.sin(a) for j in range(3)))
        normal = _norm(_cross(t, side))
        corners = []
        for sl, sw in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            corners.append(mesh.vertex(tuple(bucket['centre'][j] + t[j] * sl * bucket['half_length']
                                             + side[j] * sw * bucket['half_width'] for j in range(3)), normal))
        mesh.triangle(corners[0], corners[1], corners[2])
        mesh.triangle(corners[0], corners[2], corners[3])


def bins(bucket_list, edge):
    """Aggregate buckets into spatial bins. Each proxy keeps the site ranges of
    every bucket it swallowed, so a proxy can always name its source sites."""
    grouped = {}
    for bucket in bucket_list:
        key = tuple(int(math.floor(bucket['centre'][i] / edge)) for i in range(3))
        grouped.setdefault(key, []).append(bucket)
    out = []
    for key in sorted(grouped):
        members = grouped[key]
        # half_length + half_width bounds the bucket's oriented box on every axis,
        # so a proxy is a real bound on the sites it stands for, not an approximation.
        lo = [min(b['centre'][i] - b['half_length'] - b['half_width'] for b in members) for i in range(3)]
        hi = [max(b['centre'][i] + b['half_length'] + b['half_width'] for b in members) for i in range(3)]
        out.append({'id': 'px.%+03d%+03d%+03d' % key, 'bin': list(key),
                    'centre': [(lo[i] + hi[i]) * 0.5 for i in range(3)],
                    'half_extents': [max(1e-4, (hi[i] - lo[i]) * 0.5) for i in range(3)],
                    'sources': [{'axis': b['axis'], 'first': b['first'], 'last': b['last']} for b in members],
                    'sites': sum(b['sites'] for b in members), 'buckets': len(members)})
    return out


def proxy_box(mesh, proxy):
    """An octahedron on the proxy's extents. A deliberately plain first proxy: its
    defects are meant to be visible rather than dressed up."""
    c, e = proxy['centre'], proxy['half_extents']
    v = [mesh.vertex((c[0] + e[0], c[1], c[2]), (1, 0, 0)), mesh.vertex((c[0] - e[0], c[1], c[2]), (-1, 0, 0)),
         mesh.vertex((c[0], c[1] + e[1], c[2]), (0, 1, 0)), mesh.vertex((c[0], c[1] - e[1], c[2]), (0, -1, 0)),
         mesh.vertex((c[0], c[1], c[2] + e[2]), (0, 0, 1)), mesh.vertex((c[0], c[1], c[2] - e[2]), (0, 0, -1))]
    for a, b in ((0, 2), (2, 1), (1, 3), (3, 0)):
        mesh.triangle(v[a], v[b], v[4])
        mesh.triangle(v[b], v[a], v[5])


def compile_tier(graph, tier, scope, fixture, subject):
    if tier not in TIERS:
        raise ValueError(f'unsupported representation tier {tier!r}')
    segments = {int(k): v for k, v in fixture['wood_radial_segments'][tier].items()}
    budget = fixture['budgets'][f'{tier}_max_triangles']
    needle_low, needle_high = scope.measured('pa.needle.length_m_range')
    pitch = scope.choice('fx.needle_pitch_deg')

    wood = Mesh('wood', 'wood')
    for axis_id in sorted(graph.axes):
        axis = graph.axes[axis_id]
        if axis.order in segments:
            tube(wood, axis, segments[axis.order])
    parts, proxies, projection = [wood], [], {}

    if tier == 'near':
        sites = graph.total_site_count()
        projection = {'sites': sites, 'needle_triangles': sites * 4, 'wood_triangles': wood.triangle_count,
                      'total_triangles': sites * 4 + wood.triangle_count}
        if projection['total_triangles'] > budget:
            raise Budget(tier, projection['total_triangles'], budget, projection)
        foliage = Mesh('needles', 'foliage')
        for axis_id in sorted(graph.axes):
            for site in graph.sites(axis_id):
                needle(foliage, site)
        parts.append(foliage)
    else:
        arc = scope.choice('fx.mid_bucket_arc_m')
        bucket_list = buckets(graph, arc, needle_low, needle_high, pitch)
        if tier == 'mid':
            per = scope.choice('fx.mid_cards_per_bucket')
            projection = {'buckets': len(bucket_list), 'cluster_triangles': len(bucket_list) * per * 2,
                          'wood_triangles': wood.triangle_count,
                          'total_triangles': len(bucket_list) * per * 2 + wood.triangle_count}
            if projection['total_triangles'] > budget:
                raise Budget(tier, projection['total_triangles'], budget, projection)
            foliage = Mesh('clusters', 'foliage')
            for bucket in bucket_list:
                cards(foliage, bucket, per)
            parts.append(foliage)
            proxies = [{'id': b['id'], 'axis': b['axis'], 'first': b['first'], 'last': b['last'],
                        'sites': b['sites'], 'centre': list(b['centre']), 'tangent': list(b['tangent']),
                        'half_length': b['half_length'], 'half_width': b['half_width']} for b in bucket_list]
        else:
            proxies = bins(bucket_list, scope.choice('fx.far_bin_m'))
            projection = {'buckets': len(bucket_list), 'proxies': len(proxies),
                          'proxy_triangles': len(proxies) * 8, 'wood_triangles': wood.triangle_count,
                          'total_triangles': len(proxies) * 8 + wood.triangle_count}
            if projection['total_triangles'] > budget:
                raise Budget(tier, projection['total_triangles'], budget, projection)
            foliage = Mesh('proxies', 'foliage')
            for proxy in proxies:
                proxy_box(foliage, proxy)
            parts.append(foliage)

    return {'subject': subject, 'tier': tier, 'available': True,
            'artifact_class': graph.artifact_class, 'world_admission': False,
            'graph_digest': graph.structural_digest(), 'frame': graph.frame,
            'local_to_tree': graph.local_to_tree, 'extracted_from': graph.extracted_from,
            'parts': [p.payload() for p in parts], 'proxies': proxies,
            'counts': {'triangles': sum(p.triangle_count for p in parts),
                       'vertices': sum(p.vertex_count for p in parts),
                       'wood_triangles': wood.triangle_count,
                       'foliage_triangles': sum(p.triangle_count for p in parts) - wood.triangle_count,
                       'proxies': len(proxies),
                       'sites_represented': sum(p.get('sites', 0) for p in proxies) or graph.total_site_count(),
                       'budget': budget, 'projection': projection}}


def resolve(graph, proxy):
    """Every site a proxy stands for, resolved back out of the graph."""
    if 'sources' in proxy:
        return [sid for source in proxy['sources'] for sid in graph.site_ids(source['axis'], source['first'], source['last'])]
    return graph.site_ids(proxy['axis'], proxy['first'], proxy['last'])
