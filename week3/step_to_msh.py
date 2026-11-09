"""Tet mesh of any STEP file, checked against the CAD, written next to it as .msh. No boundary names.

Run:  .venv/bin/python week3/step_to_msh.py part.step [h_max_metres]
"""
import os
import sys

import gmsh
import numpy as np

CELLS_PER_CIRCLE = 20  # curvature sizing: on a face of radius r, h <= 2 pi r / 20


def fresh_gmsh():
    """Empty model, default options. Reuses a running gmsh: re-initialising one in a process that has
    already shown a GUI window can kill it (i.e. a notebook kernel), so never finalize between views."""
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.clear()
    gmsh.option.restoreDefaults()


def mesh_step(step, h_max):
    """Import, heal and tet-mesh a STEP file in metres, then check it. Leaves gmsh open for the caller."""
    fresh_gmsh()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setString("Geometry.OCCTargetUnit", "M")  # STEP carries its own unit (mm, inch, ...); convert it
    gmsh.model.occ.importShapes(step)
    # Fix tiny edges, sliver faces and gaps (downloaded CAD needs this). Sew loose faces into a solid only when the
    # file has no solids: sewing real solids glues touching ones into one shell (losing both) or leaves orphan faces.
    has_solids = bool(gmsh.model.occ.getEntities(3))
    gmsh.model.occ.healShapes(sewFaces=not has_solids, makeSolids=not has_solids)
    volumes = gmsh.model.occ.getEntities(3)
    assert volumes, f"{step}: no solid after healing (open shell?)"
    if len(volumes) > 1:
        gmsh.model.occ.fragment(volumes, [])  # touching solids share their faces -> one conformal mesh
    gmsh.model.occ.synchronize()
    on_solids = {t for _, t in gmsh.model.getBoundary(gmsh.model.getEntities(3), combined=False, oriented=False)}
    orphans = [t for _, t in gmsh.model.getEntities(2) if t not in on_solids]
    assert not orphans, f"faces {orphans} belong to no solid; they would be meshed as loose 2D patches"

    # Size = min(h_max, curvature size): small holes and fillets get small cells, flat regions get h_max.
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", CELLS_PER_CIRCLE)
    # Delaunay leaves the odd sliver (near-flat tet: ill-conditioned FEM matrix); Netgen's optimiser removes them.
    # linkrods.step worst minSICN 0.028 -> 0.376, for ~2x mesh time.
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.model.mesh.generate(3)

    # Checks against the exact CAD. Hard asserts catch bugs (inverted tets, lost volume); the per-face area
    # error only measures resolution, so it warns. A chord through 1/N of a circle cuts off ~(2 pi / N)^2 / 6
    # of the area; allow 3x that. Thin B-spline strips can beat this bound when h_max is coarse for them
    # (t20_data.step face 2: 11% at 5 mm, converging as h^2), hence a warning, not an abort.
    node_tags, xyz, _ = gmsh.model.mesh.getNodes()
    X = np.zeros((node_tags.max() + 1, 3))
    X[node_tags] = xyz.reshape(-1, 3)
    rtol = (2 * np.pi / CELLS_PER_CIRCLE) ** 2 / 2

    V = 0.0
    for _, vol in gmsh.model.getEntities(3):
        _, _, tet_nodes = gmsh.model.mesh.getElements(3, vol)
        T = X[tet_nodes[0].reshape(-1, 4)]
        v = np.einsum("ij,ij->i", np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), T[:, 3] - T[:, 0]) / 6
        assert (v > 0).all(), f"volume {vol}: {(v <= 0).sum()} inverted tets"
        V += v.sum()
    V_cad = sum(gmsh.model.occ.getMass(3, t) for _, t in gmsh.model.getEntities(3))
    assert abs(V - V_cad) / V_cad < rtol, f"mesh volume {V} vs CAD {V_cad}"

    face_err = {}
    for _, s in gmsh.model.getEntities(2):
        _, _, tri_nodes = gmsh.model.mesh.getElements(2, s)
        P = X[tri_nodes[0].reshape(-1, 3)]
        area = np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1).sum() / 2
        cad = gmsh.model.occ.getMass(2, s)
        face_err[s] = abs(area - cad) / cad
    worst = max(face_err.values())
    coarse = {s: e for s, e in face_err.items() if e > rtol}
    if coarse:
        print(f"WARNING: {len(coarse)} faces under-resolved (area err > {rtol:.1e}): "
              + ", ".join(f"face {s} {e:.1e}" for s, e in sorted(coarse.items(), key=lambda kv: -kv[1])[:5])
              + " -- lower h_max")

    _, tet_tags, _ = gmsh.model.mesh.getElements(3)
    quality = gmsh.model.mesh.getElementQualities(np.concatenate(tet_tags), "minSICN")
    print(f"{os.path.basename(step)}: {len(gmsh.model.getEntities(3))} volumes, {len(gmsh.model.getEntities(2))} faces, "
          f"volume {V:.4e} m^3 (CAD {V_cad:.4e}), worst face area err {worst:.1e}")
    print(f"h_max = {h_max * 1e3:g} mm: {len(node_tags)} nodes, {len(quality)} tets, "
          f"worst minSICN {quality.min():.3f} (1 = perfect, <= 0 = inverted)")


if __name__ == "__main__":
    step = sys.argv[1]
    mesh_step(step, float(sys.argv[2]) if len(sys.argv) > 2 else 0.005)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(os.path.splitext(step)[0] + ".msh")
    gmsh.finalize()
