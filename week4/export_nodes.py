"""Export mesh nodes for the neighbour-construction group: one CSV per mesh.

    python week4/export_nodes.py        (writes both CSVs and week4/nodes_summary.md)

Columns: id, x, y, z, patch, nx, ny, nz, h
  id       gmsh node tag, same as in the .msh
  patch    boundary patch name, or "interior". A node on two patches (wing root on the symmetry
           plane, box corners) gets the patch with the lowest physical tag, i.e. the wall wins.
  n        unit normal of that patch at the node, pointing out of the fluid domain (OpenFOAM
           convention: into the wing / airfoil, away from the far field); area-weighted mean of
           the adjacent boundary facets. 0, 0, 0 for interior nodes.
  h        length of the shortest mesh edge at the node (size_field.node_scale, reduce="min").
           At a boundary-layer wall node that is the first cell height.
Floats are written with 17 significant digits, so they read back bit-identical.
"""
import csv
import pathlib
import sys
from collections import Counter

import gmsh
import numpy as np
from scipy.spatial import cKDTree

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "week1"))
from mesh_check import DIM, FACES, extract_gmsh  # noqa: E402
from size_field import _local, node_scale  # noqa: E402

CASES = [  # (.msh, CSV). 2-D: the gmsh mesh, not the extruded OpenFOAM mesh (every node twice)
    (HERE / "onera_m6_wing.msh", HERE / "onera_m6_nodes.csv"),
    (HERE / "2d" / "onera_d_airfoil_2d.msh", HERE / "2d" / "onera_d_2d_nodes.csv"),
]
DUPLICATE_TOL = 1e-12  # m


def read(msh):
    gmsh.initialize(["-v", "0"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(msh))
        return extract_gmsh()
    finally:
        gmsh.finalize()


def facet_normals(xyz, facets):
    """(F, k) local indices -> (F, 3) facet normals, length = facet length (2-D) or area (3-D)."""
    p = xyz[facets]
    if facets.shape[1] == 2:  # 2-D edge in the xy plane
        t = p[:, 1] - p[:, 0]
        return np.stack([t[:, 1], -t[:, 0], np.zeros(len(t))], axis=1)
    if facets.shape[1] == 3:
        return 0.5 * np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    raise ValueError(f"boundary facets with {facets.shape[1]} nodes not handled")


def patches_and_normals(node_tags, xyz, elements, boundary):
    """Per node: patch name ("interior" if none) and outward unit normal (0 if interior)."""
    n = len(node_tags)
    patch = np.full(n, "interior", dtype=object)
    normal = np.zeros((n, 3))
    on_patches = np.zeros(n, int)
    boundary = {name: _local(node_tags, f) for name, f in boundary.items()}

    # cell centroid behind each boundary facet, to orient the facet normal outward
    wanted = {k for f in boundary.values() for k in map(tuple, np.sort(f, axis=1))}
    behind = {}
    top_dim = max(DIM[t] for t in elements)
    for elem_type, conn in elements.items():
        if DIM[elem_type] != top_dim:
            continue
        cells = _local(node_tags, conn)
        centroid = xyz[cells].mean(axis=1)
        for facet in FACES[elem_type]:
            for i, key in enumerate(map(tuple, np.sort(cells[:, facet], axis=1))):
                if key in wanted:
                    behind[key] = centroid[i]

    for name, facets in reversed(boundary.items()):  # reversed: lowest physical tag written last, wins
        fn = facet_normals(xyz, facets)
        inward = np.array([behind[k] for k in map(tuple, np.sort(facets, axis=1))])
        fn *= np.sign(((xyz[facets].mean(axis=1) - inward) * fn).sum(axis=1))[:, None]
        acc = np.zeros((n, 3))
        for j in range(facets.shape[1]):
            np.add.at(acc, facets[:, j], fn)
        nodes = np.unique(facets)
        patch[nodes] = name
        normal[nodes] = acc[nodes] / np.linalg.norm(acc[nodes], axis=1, keepdims=True)
        on_patches[nodes] += 1
    return patch, normal, int((on_patches > 1).sum())


def export(msh, out):
    m = read(msh)
    tags, xyz = m["node_tags"].astype(np.int64), m["coords"].reshape(-1, 3)
    h = node_scale(tags, xyz, m["elements"], reduce="min")
    patch, normal, shared = patches_and_normals(tags, xyz, m["elements"], m["boundary"])
    duplicates = len(cKDTree(xyz).query_pairs(DUPLICATE_TOL))

    assert not np.isnan(h).any(), "node without a domain element"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "x", "y", "z", "patch", "nx", "ny", "nz", "h"])
        g = "%.17g".__mod__
        for t, p, pa, nv, hi in zip(tags, xyz, patch, normal, h):
            w.writerow([int(t), g(p[0]), g(p[1]), g(p[2]), pa, g(nv[0]), g(nv[1]), g(nv[2]), g(hi)])

    # read back: 17 digits must reproduce every float exactly
    back = np.loadtxt(out, delimiter=",", skiprows=1, usecols=(1, 2, 3, 5, 6, 7, 8))
    assert np.array_equal(back[:, :3], xyz) and np.array_equal(back[:, 6], h), "CSV round-trip lost digits"
    print(f"{out.relative_to(HERE.parent)}: {len(tags)} nodes, "
          f"{duplicates} duplicate pairs (distance < {DUPLICATE_TOL:g} m)")
    return dict(out=out, n=len(tags), duplicates=duplicates, shared=shared,
                patches=Counter(patch), patch=patch, h=h, z_max=np.abs(xyz[:, 2]).max())


def summary(results):
    lines = ["# Node export summary", "",
             "Generated by `python week4/export_nodes.py`. Columns and conventions: see the script docstring.", ""]
    for r in results:
        h = r["h"]
        lines += [f"## `{r['out'].relative_to(HERE)}`", "",
                  "| Patch | Nodes | h min |", "|---|---|---|"]
        lines += [f"| {p} | {c:,} | {h[r['patch'] == p].min():.4g} |" for p, c in sorted(r["patches"].items())]
        lines += [f"| **total** | **{r['n']:,}** | |", "",
                  "| h min | h median | h max |", "|---|---|---|",
                  f"| {h.min():.4g} | {np.median(h):.4g} | {h.max():.4g} |", "",
                  f"- Duplicate / near-duplicate node pairs (distance < {DUPLICATE_TOL:g} m): {r['duplicates']}",
                  f"- Nodes on two or more patches (assigned to the lowest physical tag): {r['shared']:,}",
                  f"- max abs(z): {r['z_max']:g}", ""]
    (HERE / "nodes_summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    summary([export(msh, out) for msh, out in CASES])
