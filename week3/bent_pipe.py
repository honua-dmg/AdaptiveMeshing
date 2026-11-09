"""90-degree bent pipe with straight legs, written to pipe.step.

Run:  /Applications/FreeCAD.app/Contents/Resources/bin/freecadcmd week3/bent_pipe.py

FreeCAD works in mm, so the STEP file is in mm; scale to metres on gmsh import.

Layout (z = 0 plane, bend centre at the origin):

        end_b
   (-L, R) ===========  (0, R)
                         \\  bend: centreline radius R, 90 degrees
                          ||  (R, 0)
                          ||
                          ||  leg a, along +y
                          ||
                       (R, -L)
                        end_a
"""
import math
import os

import Part
from FreeCAD import Vector as V

# Every dimension is a parameter (mm).
R_OUT = 50.0    # outer pipe radius
R_IN = 40.0     # inner pipe radius -> wall thickness 10
R_BEND = 150.0  # bend radius of the centreline (1.5 x nominal diameter, a "long radius" elbow)
L = 300.0       # straight leg length; long enough for the Lame / log-profile checks to hold mid-leg


def tube(r):
    """Solid of radius r swept along the centreline: leg a + 90 degree bend + leg b."""
    leg_a = Part.makeCylinder(r, L, V(R_BEND, -L, 0), V(0, 1, 0))
    bend = Part.makeTorus(R_BEND, r, V(0, 0, 0), V(0, 0, 1), -180, 180, 90)
    leg_b = Part.makeCylinder(r, L, V(0, R_BEND, 0), V(-1, 0, 0))
    return leg_a.fuse([bend, leg_b])


pipe = tube(R_OUT).cut(tube(R_IN))

# Pappus: annulus area x centreline length.
exact = math.pi * (R_OUT**2 - R_IN**2) * (2 * L + math.pi * R_BEND / 2)
assert len(pipe.Solids) == 1, f"expected 1 solid, got {len(pipe.Solids)}"
assert abs(pipe.Volume - exact) / exact < 1e-6, f"volume {pipe.Volume} != {exact}"

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipe.step")
pipe.exportStep(out)
print(f"wrote {out}: volume {pipe.Volume:.1f} mm^3 (exact {exact:.1f}), {len(pipe.Faces)} faces")
