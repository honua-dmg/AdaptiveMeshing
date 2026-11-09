"""Pipe mesh for the FEM solver: step_to_msh's generic mesh + named boundary groups, written to pipe.msh.

Run:  .venv/bin/python week3/mesh_pipe.py [h_metres]

Mesh contract (what the FEM solver reads):
  units metres, gmsh 4.1 format, 4-node tets
  volume group   wall
  surface groups inner, outer, end_a, end_b
"""
import os
import sys

import gmsh
import numpy as np

from step_to_msh import mesh_step

H = float(sys.argv[1]) if len(sys.argv) > 1 else 0.005  # 5 mm = 2 cells through the 10 mm wall
HERE = os.path.dirname(os.path.abspath(__file__))

mesh_step(os.path.join(HERE, "pipe.step"), H)


# Tag faces by shape, not by tag number (numbers shift whenever the CAD changes).
# Largest principal curvature: 0 on a flat end, 1/r on a cylinder or torus wall,
# so the inner wall (smaller r) curves more than the outer one.
def max_curvature(tag):
    lo, hi = gmsh.model.getParametrizationBounds(2, tag)
    kmax, kmin, _, _ = gmsh.model.getPrincipalCurvatures(tag, (np.asarray(lo) + np.asarray(hi)) / 2)
    return max(abs(kmax[0]), abs(kmin[0]))


surfaces = [t for _, t in gmsh.model.getEntities(2)]
k = {t: max_curvature(t) for t in surfaces}
flat = [t for t in surfaces if k[t] < 1e-6]
curved = [t for t in surfaces if t not in flat]
k_in = max(k[t] for t in curved)
inner = [t for t in curved if np.isclose(k[t], k_in, rtol=1e-3)]
outer = [t for t in curved if t not in inner]
assert len(flat) == 2, f"expected 2 flat end faces, got {len(flat)}"
assert np.allclose([k[t] for t in outer], k[outer[0]], rtol=1e-3), "more than two wall radii"
# end_a closes leg a (lowest y); end_b closes leg b.
end_a, end_b = sorted(flat, key=lambda t: gmsh.model.occ.getCenterOfMass(2, t)[1])
print(f"radii from curvature: inner {1 / k_in * 1e3:.1f} mm, outer {1 / k[outer[0]] * 1e3:.1f} mm")

for name, tags in {"inner": inner, "outer": outer, "end_a": [end_a], "end_b": [end_b]}.items():
    gmsh.model.addPhysicalGroup(2, tags, name=name)
(_, vol), = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, [vol], name="wall")

gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.write(os.path.join(HERE, "pipe.msh"))
gmsh.finalize()
