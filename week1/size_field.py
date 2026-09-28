"""Prescribed target size field -> gmsh -> realized-vs-requested diagnostic.

Requested resolution is h_test(x), an analytic field we choose. Realized resolution is
whatever gmsh actually built. diagnose() reports the multiplicative ratio and its log:

    h_mesh_i = mean length of the domain edges incident to node i
    r_i      = h_mesh_i / (h_target_i + eps)
    delta_i  = log(r_i)                        delta ~ 0  <=>  mesh matches request

Two ways to hand h_test to gmsh, so a mismatch can be attributed:
  "matheval" - closed form as a gmsh expression. No sampling error. Ground truth.
  "postview" - h_test sampled on a background mesh, consumed as view data. Same
               mechanism a PDE error estimator will use later.
"""
import pathlib
import sys

import numpy as np

from mesh_check import DIM, FACES, extract_gmsh

# gmsh takes min() over every size source then clamps to [MeshSizeMin, MeshSizeMax], so a
# nonzero floor or a leftover point size silently overrides the field. Kill all of them.
FIELD_ONLY = {
    "Mesh.MeshSizeFromPoints": 0,
    "Mesh.MeshSizeFromCurvature": 0,
    "Mesh.MeshSizeExtendFromBoundary": 0,
    "Mesh.MeshSizeMin": 0,
    "Mesh.MeshSizeMax": 1e22,
}


def h_test(coords, hmin=0.01, hmax=0.1, xc=(0.5, 0.5, 0.0), sigma=0.15):
    """hmin + (hmax - hmin) * (1 - exp(-||x - xc||^2 / (2 sigma^2))); (N,3) -> (N,)."""
    d2 = ((np.asarray(coords, float).reshape(-1, 3) - np.asarray(xc, float)) ** 2).sum(1)
    return hmin + (hmax - hmin) * (1.0 - np.exp(-d2 / (2.0 * sigma ** 2)))


def _matheval_expr(hmin, hmax, xc, sigma):
    """h_test as a gmsh MathEval string. Must stay in sync with h_test above."""
    d2 = f"((x-({xc[0]!r}))^2+(y-({xc[1]!r}))^2+(z-({xc[2]!r}))^2)"
    return f"({hmin!r})+({(hmax - hmin)!r})*(1-Exp(-{d2}/(2*({sigma!r})^2)))"


def _edges(elem_type):
    """Local node index pairs of an element's edges, derived from its facets.

    A facet list gives edges directly in 2D (facets *are* edges) and, in 3D, every edge
    appears as a consecutive pair around some face. Saves keeping a second table.
    """
    return {tuple(sorted((f[i], f[(i + 1) % len(f)])))
            for f in FACES[elem_type] for i in range(len(f))}


def _unit_square(dim):
    """Unit square in z=0 (dim 2) or unit cube (dim 3). Point sizes left at 0: field only."""
    import gmsh

    if dim == 2:
        gmsh.model.occ.add_rectangle(0, 0, 0, 1, 1)
    else:
        gmsh.model.occ.add_box(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()


def _mesh_in_session(mode="matheval", dim=2, bg_size=None, name="target", **params):
    """Geometry + size field + generate, inside an already-open gmsh session.

    Leaves the new model current, so the caller can gmsh.write() it or read views off it.
    build() wraps this with init/finalize; show() keeps the session open for the GUI.
    """
    import gmsh

    view = None
    if mode == "postview":
        gmsh.model.add(f"{name}__bg")
        _unit_square(dim)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", bg_size or params.get("hmin", 0.01))
        gmsh.model.mesh.generate(dim)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        view = gmsh.view.add(f"{name}__h_target_bg")
        gmsh.view.addHomogeneousModelData(
            view, 0, f"{name}__bg", "NodeData", tags, h_test(coords, **params))

    gmsh.model.add(name)
    _unit_square(dim)

    if mode == "matheval":
        f = gmsh.model.mesh.field.add("MathEval")
        gmsh.model.mesh.field.setString(f, "F", _matheval_expr(
            params.get("hmin", 0.01), params.get("hmax", 0.1),
            params.get("xc", (0.5, 0.5, 0.0)), params.get("sigma", 0.15)))
    elif mode == "postview":
        f = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(f, "ViewTag", view)
    else:
        raise ValueError(f"mode must be 'matheval' or 'postview', got {mode!r}")

    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    for option, value in FIELD_ONLY.items():
        gmsh.option.setNumber(option, value)

    gmsh.model.mesh.generate(dim)
    return view


def build(mode="matheval", dim=2, bg_size=None, **params):
    """Mesh the unit domain under h_test(**params). Returns extract_gmsh()'s dict.

    mode "postview" first meshes a throwaway background model at bg_size, samples h_test
    on its nodes, and feeds that as view data -- the shape of a real adaptive step. bg_size
    defaults to hmin so the grid resolves the whole requested range; a coarser grid cannot
    represent the well and biases delta upward, which is sampling error, not gmsh.
    """
    import gmsh

    gmsh.initialize(["-noenv"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        _mesh_in_session(mode=mode, dim=dim, bg_size=bg_size, **params)
        return extract_gmsh()
    finally:
        gmsh.finalize()


def _local(node_tags, conn):
    """Element connectivity given in node tags -> row indices into the coords array."""
    node_tags = np.asarray(node_tags)
    order = np.argsort(node_tags)
    return order[np.searchsorted(node_tags[order], np.asarray(conn))]


def node_scale(node_tags, coords, elements, reduce="mean"):
    """Local mesh scale per node from its incident domain edges; nan where it has none.

    reduce="mean" averages them; reduce="min" takes the shortest. Measured, both are
    unbiased on a uniform mesh (median delta -0.000 and -0.001), so neither carries a
    constant offset. They diverge only where the target field is graded, and the size of
    the gap measures how spread a node's incident edge lengths are. Run both: agreement
    means the local edges are uniform, divergence localises the transition regions.
    """
    node_tags = np.asarray(node_tags)
    coords = np.asarray(coords, float).reshape(-1, 3)

    if reduce not in ("mean", "min"):
        raise ValueError(f"reduce must be 'mean' or 'min', got {reduce!r}")
    accumulate = np.add.at if reduce == "mean" else np.fmin.at
    n = len(node_tags)
    total = np.zeros(n) if reduce == "mean" else np.full(n, np.inf)
    count = np.zeros(n)

    top_dim = max(DIM[t] for t in elements)
    for elem_type, conn in elements.items():
        if DIM[elem_type] != top_dim:
            continue  # lower-dim elements are boundary copies; their edges double-count
        local = _local(node_tags, conn)
        for a, b in _edges(elem_type):
            length = np.linalg.norm(coords[local[:, a]] - coords[local[:, b]], axis=1)
            for side in (local[:, a], local[:, b]):
                accumulate(total, side, length)
                np.add.at(count, side, 1.0)

    if reduce == "min":
        return np.where(count > 0, total, np.nan)
    return np.divide(total, count, out=np.full(n, np.nan), where=count > 0)


def diagnose(node_tags, coords, elements, boundary=None, eps=1e-12, bins=6,
             reduce="mean", **params):
    """Per-node realized-vs-requested report. **params go to h_test.

    Returns {"delta", "r", "h_mesh", "h_target", "radius", summary scalars, "profile"}.
    eps guards h_target -> 0; inert for this field (h_target >= hmin) but AI-generated
    targets can request zero.
    """
    coords = np.asarray(coords, float).reshape(-1, 3)
    h_mesh = node_scale(node_tags, coords, elements, reduce=reduce)
    h_target = h_test(coords, **params)
    r = h_mesh / (h_target + eps)
    delta = np.log(r)

    xc = np.asarray(params.get("xc", (0.5, 0.5, 0.0)), float)
    radius = np.linalg.norm(coords - xc, axis=1)
    ok = np.isfinite(delta)

    # trend: is realized size actually growing with distance from xc, as requested?
    edges = np.quantile(radius[ok], np.linspace(0, 1, bins + 1))
    which = np.clip(np.searchsorted(edges, radius, side="right") - 1, 0, bins - 1)
    profile = [dict(r=float(np.median(radius[m])), n=int(m.sum()),
                    h_target=float(np.median(h_target[m])), h_mesh=float(np.median(h_mesh[m])),
                    delta=float(np.median(delta[m])))
               for b in range(bins) if (m := ok & (which == b)).any()]

    q1, q3 = np.quantile(delta[ok], [0.25, 0.75])
    return dict(delta=delta, r=r, h_mesh=h_mesh, h_target=h_target, radius=radius,
                nodes=int(ok.sum()), unscaled=int((~ok).sum()),
                median_delta=float(np.median(delta[ok])), iqr_delta=float(q3 - q1),
                mean_abs_delta=float(np.abs(delta[ok]).mean()),
                median_r=float(np.median(r[ok])), profile=profile)


LIST_TYPE = {2: "ST", 3: "SQ", 4: "SS", 5: "SH"}  # gmsh list-data names, scalar per node


def add_view(name, node_tags, coords, elements, values):
    """Standalone post-processing view of a nodal field, drawn on the domain elements.

    ListData carries its own coordinates, so the view does not depend on the model it came
    from still being current -- which is what lets one GUI window hold every case at once.
    Gmsh renders only the current model, so views are the way to compare meshes side by side.
    """
    import gmsh

    coords = np.asarray(coords, float).reshape(-1, 3)
    values = np.asarray(values, float)
    ok = np.isfinite(values)
    values = np.where(ok, values, np.median(values[ok]) if ok.any() else 0.0)

    tag = gmsh.view.add(name)
    top_dim = max(DIM[t] for t in elements)
    for elem_type, conn in elements.items():
        if DIM[elem_type] != top_dim:
            continue
        local = _local(node_tags, conn)
        pts = coords[local]                                  # (M, nodes per elem, 3)
        rows = np.concatenate([pts.transpose(0, 2, 1).reshape(len(pts), -1),
                               values[local]], axis=1)       # x1..xk, y1..yk, z1..zk, v1..vk
        gmsh.view.addListData(tag, LIST_TYPE[elem_type], len(rows), rows.ravel())
    return tag


def show(cases=None, modes=("matheval", "postview"), dim=2,
         fields=("h_target", "h_mesh", "delta"), outdir="meshes", run_gui=True):
    """Mesh every case, save each .msh, open one gmsh GUI holding all of them as views.

    All views start hidden except the first; toggle them in the GUI's view list to flip
    between cases without regenerating anything. In 3D the tets fill the volume, so use the
    GUI's clipping planes (Tools > Clipping) to see inside. Returns [(label, case, report)]
    so the numbers line up with what is on screen.
    """
    import gmsh

    cases = CASES if cases is None else cases
    directory = pathlib.Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)

    gmsh.initialize(["-noenv"])
    rows = []
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.SaveAll", 1)   # no physical groups here; save anyway
        for i, case in enumerate(cases):
            for mode in modes:
                label = f"c{i}_{mode}"
                _mesh_in_session(mode=mode, dim=dim, name=label, **case)
                gmsh.write(str(directory / f"{label}.msh"))

                mesh = extract_gmsh()
                d = diagnose(**mesh, **case)
                tail = (f"hmin {case['hmin']:g}  hmax {case['hmax']:g}  "
                        f"sigma {case['sigma']:g}")
                for field in fields:
                    add_view(f"{label}  {field}   {tail}",
                             mesh["node_tags"], mesh["coords"], mesh["elements"], d[field])
                rows.append((label, case, d))

        # The last model meshed stays current, and gmsh draws its mesh on top of every view
        # -- an unrelated case overlaid on the one you are looking at, which reads as "the
        # mesh never changed". Views carry their own wireframe via ShowElement, so turn all
        # model and geometry drawing off and let the views be the only thing on screen.
        for option in ("Mesh.SurfaceEdges", "Mesh.SurfaceFaces", "Mesh.VolumeEdges",
                       "Mesh.VolumeFaces", "Mesh.Points", "Mesh.Lines",
                       "Geometry.Points", "Geometry.Curves", "Geometry.Surfaces"):
            gmsh.option.setNumber(option, 0)

        for tag in gmsh.view.getTags():
            index = gmsh.view.getIndex(tag)
            gmsh.option.setNumber(f"View[{index}].ShowElement", 1)   # draw the mesh edges
            gmsh.option.setNumber(f"View[{index}].Visible", 1 if index == 0 else 0)

        if run_gui and "close" not in sys.argv:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()
    return rows


CASES = [
    dict(hmin=0.01, hmax=0.10, sigma=0.15, xc=(0.5, 0.5, 0.0)),   # baseline well at centre
    dict(hmin=0.02, hmax=0.10, sigma=0.05, xc=(0.5, 0.5, 0.0)),   # tight well: steep gradient
    dict(hmin=0.01, hmax=0.20, sigma=0.30, xc=(0.5, 0.5, 0.0)),   # wide well, 20x size range
    dict(hmin=0.02, hmax=0.08, sigma=0.12, xc=(0.0, 0.0, 0.0)),   # well on a corner
    dict(hmin=0.05, hmax=0.05, sigma=0.15, xc=(0.5, 0.5, 0.0)),   # uniform control: no trend
]


def sweep(cases=CASES, modes=("matheval", "postview"), dim=2):
    """Run every case through every mode; return [(mode, case, report)]."""
    return [(mode, case, diagnose(**build(mode=mode, dim=dim, **case), **case))
            for case in cases for mode in modes]


def report(rows):
    head = f"{'mode':9} {'hmin':>5} {'hmax':>5} {'sigma':>5} {'nodes':>6} " \
           f"{'med delta':>9} {'IQR':>6} {'mean|d|':>7} {'med r':>6}"
    print(head, "-" * len(head), sep="\n")
    for mode, case, d in rows:
        print(f"{mode:9} {case['hmin']:5.3f} {case['hmax']:5.3f} {case['sigma']:5.3f} "
              f"{d['nodes']:6d} {d['median_delta']:9.3f} {d['iqr_delta']:6.3f} "
              f"{d['mean_abs_delta']:7.3f} {d['median_r']:6.3f}")


def _selfcheck():
    # h_test: hmin at the centre, -> hmax far away, monotone in radius
    assert np.isclose(h_test([[0.5, 0.5, 0]])[0], 0.01)
    assert np.isclose(h_test([[0.5, 0.5, 0]], hmin=0.2, hmax=0.2)[0], 0.2)
    far = h_test([[0.5, 0.5, 0], [0.6, 0.5, 0], [1.5, 0.5, 0]])
    assert (np.diff(far) > 0).all() and far[-1] < 0.1

    # edges derived from facets: 3 for a tri, 4 for a quad, 6 for a tet, 12 for a hex
    assert {len(_edges(t)) for t in (2,)} == {3} and len(_edges(3)) == 4
    assert len(_edges(4)) == 6 and len(_edges(5)) == 12

    # node_scale on a unit right triangle: legs 1, hypotenuse sqrt(2)
    tags, xy, el = [1, 2, 3], [[0, 0, 0], [1, 0, 0], [0, 1, 0]], {2: np.array([[1, 2, 3]])}
    scale = node_scale(tags, xy, el)
    assert np.allclose(scale[0], 1.0) and np.allclose(scale[1:], (1 + 2 ** 0.5) / 2), scale
    shortest = node_scale(tags, xy, el, reduce="min")
    assert np.allclose(shortest, 1.0), shortest       # every node touches a leg of length 1
    assert (node_scale(tags, xy, el, reduce="min")
            <= node_scale(tags, xy, el)).all()        # min can never exceed mean

    # add_view: the ListData payload lands, one entry per domain element
    import gmsh
    uniform_case = dict(hmin=0.05, hmax=0.05, sigma=0.15)
    mesh = build(**uniform_case)
    gmsh.initialize(["-noenv"])
    try:
        d = diagnose(**mesh, **uniform_case)
        tag = add_view("t", mesh["node_tags"], mesh["coords"], mesh["elements"], d["delta"])
        kind, numele, data = gmsh.view.getListData(tag)
        n_tri = len(mesh["elements"][2])
        assert list(kind) == ["ST"] and list(numele) == [n_tri], (kind, numele, n_tri)
        assert len(data[0]) == n_tri * 12, len(data[0])   # 9 coords + 3 values per triangle
    finally:
        gmsh.finalize()

    # show(): every case in one session, a .msh per case on disk, no GUI under test
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        shown = show(cases=CASES[:2], outdir=tmp, run_gui=False)
        assert [label for label, _, _ in shown] == [
            "c0_matheval", "c0_postview", "c1_matheval", "c1_postview"], shown
        written = sorted(f.name for f in pathlib.Path(tmp).glob("*.msh"))
        assert len(written) == 4 and all(f.stat().st_size > 0
                                         for f in pathlib.Path(tmp).glob("*.msh")), written

    # calibration that makes every other delta readable: on a uniform target neither
    # estimator may show a median offset, so a non-zero delta elsewhere is the mesh, not us
    flat = dict(hmin=0.05, hmax=0.05, sigma=0.15)
    flat_mesh = build(**flat)
    for reduce in ("mean", "min"):
        offset = diagnose(**flat_mesh, reduce=reduce, **flat)["median_delta"]
        assert abs(offset) < 0.02, (reduce, offset)

    rows = sweep()
    by_mode = {}
    for mode, case, d in rows:
        assert d["nodes"] > 0 and d["unscaled"] == 0, (mode, case, d["nodes"], d["unscaled"])
        assert np.isfinite(d["median_delta"]) and np.isfinite(d["iqr_delta"])
        assert len(d["profile"]) > 1
        by_mode.setdefault(mode, []).append((case, d))

    for mode, results in by_mode.items():
        uniform = [d for case, d in results if case["hmin"] == case["hmax"]][0]
        graded = [d for case, d in results if case["sigma"] == 0.15 and case["hmin"] == 0.01][0]
        # uniform target -> uniform mesh -> delta nearly constant; graded target spreads more
        assert uniform["iqr_delta"] < 0.1, (mode, uniform["iqr_delta"])
        assert graded["iqr_delta"] > uniform["iqr_delta"], mode
        # requested trend realized: median h_mesh grows from the innermost to outermost bin
        inner, outer = graded["profile"][0], graded["profile"][-1]
        assert outer["h_mesh"] > 2 * inner["h_mesh"], (mode, inner, outer)

    report(rows)
    print("\nradius profile, matheval baseline case:")
    baseline = next(d for mode, case, d in rows if mode == "matheval" and case["sigma"] == 0.15
                    and case["hmin"] == 0.01)
    print(f"{'r':>6} {'n':>6} {'h_target':>9} {'h_mesh':>8} {'delta':>7}")
    for row in baseline["profile"]:
        print(f"{row['r']:6.3f} {row['n']:6d} {row['h_target']:9.4f} "
              f"{row['h_mesh']:8.4f} {row['delta']:7.3f}")
    print("\nselfcheck ok")


if __name__ == "__main__":
    _selfcheck()
