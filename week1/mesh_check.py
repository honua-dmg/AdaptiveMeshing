"""Pass/fail sanity report for extracted mesh data.

check_mesh() takes plain arrays (gmsh node/element tags), extract_gmsh() builds them
from the current gmsh model. Linear elements only: line/tri/quad/tet/hex/point.
"""
from collections import Counter, defaultdict

import numpy as np

# gmsh element type -> dimension
DIM = {15: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 3}

# gmsh element type -> local node indices of each facet (edges in 2D, faces in 3D)
FACES = {
    2: [(0, 1), (1, 2), (2, 0)],
    3: [(0, 1), (1, 2), (2, 3), (3, 0)],
    4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    5: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
}


def _check(ok, **info):
    """One entry of the report: pass flag plus whatever numbers explain it."""
    return {"pass": bool(ok), **info}


def _facet_keys(rows):
    """Facet ids that ignore node order: each row's node tags, sorted, as a tuple."""
    return [tuple(r) for r in np.sort(np.asarray(rows), axis=1)]


def _count_repeated_nodes(conn):
    """Elements using the same node twice (zero area/volume)."""
    neighbour_gap = np.diff(np.sort(conn, axis=1), axis=1)
    return int((neighbour_gap == 0).any(axis=1).sum())


def _boundary_check(domain, boundary):
    """Every tagged facet must be a facet of some domain element."""
    facets = Counter()  # facet key -> number of domain elements sharing it
    for elem_type, conn in domain.items():
        for facet in FACES.get(elem_type, []):
            facets.update(_facet_keys(conn[:, facet]))

    tagged, per_tag = set(), {}
    for tag, facet_nodes in boundary.items():
        keys = _facet_keys(facet_nodes)
        unknown = [k for k in keys if k not in facets]
        tagged.update(keys)
        per_tag[tag] = _check(keys and not unknown, facets=len(keys), not_a_domain_face=len(unknown))

    outer = {f for f, shared_by in facets.items() if shared_by == 1}
    return _check(
        all(t["pass"] for t in per_tag.values()), tags=per_tag,
        untagged_outer_facets=len(outer - tagged))  # info only: partial tagging is legal


def check_mesh(node_tags, coords, elements, boundary=None, expected_bounds=None, tol=1e-9):
    """node_tags (N,), coords (N,3), elements {gmsh_type: (M,k) node tags},
    boundary {tag: (F,k) facet node tags, one facet type per tag} or None to skip,
    expected_bounds ((xmin,ymin,zmin), (xmax,ymax,zmax)) or None.
    Returns {"pass": bool, "checks": {name: {"pass": bool, ...}}}."""
    node_tags = np.asarray(node_tags)
    coords = np.asarray(coords, float).reshape(-1, 3)
    elements = {t: np.asarray(c) for t, c in elements.items()}
    unsupported = set(elements) - set(DIM)
    if unsupported:
        raise ValueError(f"unsupported element types: {sorted(unsupported)}")

    # domain = elements of the highest dimension present; lower-dim ones are boundary/embedded
    top_dim = max((DIM[t] for t in elements), default=-1)
    domain = {t: c for t, c in elements.items() if DIM[t] == top_dim}
    checks = {}

    n_nodes, n_unique = len(node_tags), len(np.unique(node_tags))
    checks["node_count"] = _check(
        n_nodes > 0 and n_nodes == len(coords) and n_unique == n_nodes,
        count=n_nodes, coord_rows=len(coords), duplicate_tags=n_nodes - n_unique)

    n_domain = sum(len(c) for c in domain.values())
    checks["element_count"] = _check(
        n_domain > 0, total=sum(len(c) for c in elements.values()), domain=n_domain,
        by_type={t: len(c) for t, c in elements.items()})

    lo = coords.min(axis=0) if n_nodes else np.full(3, np.nan)
    hi = coords.max(axis=0) if n_nodes else np.full(3, np.nan)
    as_expected = expected_bounds is None or (
        np.allclose(lo, expected_bounds[0], atol=tol) and np.allclose(hi, expected_bounds[1], atol=tol))
    checks["coordinate_bounds"] = _check(
        n_nodes > 0 and np.isfinite(coords).all() and as_expected, min=lo.tolist(), max=hi.tolist())

    used = np.concatenate([c.ravel() for c in elements.values()]) if elements else np.array([], int)
    missing = np.setdiff1d(used, node_tags)   # elements pointing at nodes that don't exist
    orphan = np.setdiff1d(node_tags, used)    # nodes no element uses
    degenerate = sum(_count_repeated_nodes(c) for c in elements.values())
    checks["connectivity"] = _check(
        not (missing.size or orphan.size or degenerate),
        missing_nodes=missing.size, missing_sample=missing[:10].tolist(),
        orphan_nodes=orphan.size, orphan_sample=orphan[:10].tolist(), degenerate_elements=degenerate)

    # ponytail: facet keys counted in python dicts; vectorise per facet width at 1e7+ elements
    checks["boundary_tags"] = _check(True, skipped=True) if boundary is None \
        else _boundary_check(domain, boundary)

    return {"pass": all(c["pass"] for c in checks.values()), "checks": checks}


def extract_gmsh():
    """Build check_mesh() args from the current gmsh model (call after generate() or open())."""
    import gmsh

    def connectivity(dim=-1, tag=-1):
        """{gmsh element type: (M, nodes per element) array of node tags}"""
        types, _, node_lists = gmsh.model.mesh.getElements(dim, tag)
        return {int(t): nodes.reshape(-1, gmsh.model.mesh.getElementProperties(t)[3])  # [3] = nodes/elem
                for t, nodes in zip(types, node_lists)}

    node_tags, coords, _ = gmsh.model.mesh.getNodes()

    # physical groups one dimension below the model are the boundary tags
    by_name_type = defaultdict(list)
    for dim, group in gmsh.model.getPhysicalGroups(gmsh.model.getDimension() - 1):
        name = gmsh.model.getPhysicalName(dim, group) or group
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, group):
            for elem_type, facets in connectivity(dim, entity).items():
                by_name_type[name, elem_type].append(facets)

    # a group holding two facet types (mixed tri/quad surface) splits into "name:type"
    shared = {n for n, count in Counter(n for n, _ in by_name_type).items() if count > 1}
    boundary = {(f"{n}:{t}" if n in shared else n): np.concatenate(arrays)
                for (n, t), arrays in by_name_type.items()}

    return dict(node_tags=node_tags, coords=coords, elements=connectivity(), boundary=boundary or None)


def _selfcheck():
    n, xy = [1, 2, 3, 4], [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    el = {2: [[1, 2, 3], [1, 3, 4]]}
    bd = {"bottom": [[1, 2]], "right": [[2, 3]], "top": [[3, 4]], "left": [[4, 1]]}
    ok = check_mesh(n, xy, el, bd, expected_bounds=([0, 0, 0], [1, 1, 0]))
    assert ok["pass"] and ok["checks"]["boundary_tags"]["untagged_outer_facets"] == 0, ok
    assert check_mesh(n, xy, el)["checks"]["boundary_tags"]["skipped"]
    part = check_mesh(n, xy, el, {k: v for k, v in bd.items() if k != "left"})
    assert part["pass"] and part["checks"]["boundary_tags"]["untagged_outer_facets"] == 1
    assert not check_mesh(n, xy, el, {"x": [[2, 4]]})["checks"]["boundary_tags"]["pass"]  # not a face
    assert not check_mesh(n, xy, el, {"x": np.empty((0, 2), int)})["pass"]                # empty tag
    assert not check_mesh(n, xy, {2: [[1, 2, 9], [1, 3, 4]]})["checks"]["connectivity"]["pass"]  # dangling ref
    assert not check_mesh(n + [5], xy + [[5, 5, 0]], el)["checks"]["connectivity"]["pass"]       # orphan node
    assert not check_mesh(n, xy, {2: [[1, 2, 2], [1, 3, 4]]})["checks"]["connectivity"]["pass"] # degenerate
    assert not check_mesh(n, xy[:3] + [[np.nan, 0, 0]], el)["checks"]["coordinate_bounds"]["pass"]
    assert not check_mesh(n, xy, el, expected_bounds=([0, 0, 0], [2, 1, 0]))["pass"]
    assert not check_mesh([1, 1, 3, 4], xy, el)["checks"]["node_count"]["pass"]              # dup tags
    print("selfcheck ok")


if __name__ == "__main__":
    _selfcheck()
