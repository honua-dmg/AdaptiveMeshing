# ONERA M6 half-wing: volume mesh

**281,604 tetrahedral cells, 57,195 points. checkMesh: `Mesh OK`.**

Unstructured tet mesh of the fluid around the half-wing (box minus wing, symmetry plane at y = 0), generated from the
supplied `onera_m6_wing.geo` unchanged. Items marked ★ need review before solving.

## Files

| File | What |
|---|---|
| `onera_m6_wing.msh` | gmsh mesh, format 2.2 ASCII, 16 MB |
| `case/constant/polyMesh/` | OpenFOAM mesh (from `gmshToFoam`) |
| `gmsh.log`, `case/log.gmshToFoam`, `case/log.checkMesh` | full tool output |

| Patch | Faces |
|---|---|
| `wing` | 38,902 |
| `symmetry` | 6,034 |
| `farfield` | 804 |
| `inlet` | 122 |
| `outlet` | 122 |

## Table 1: Topology (checkMesh)

| Check | Result |
|---|---|
| Boundary definition | OK |
| Cell-to-face addressing | OK |
| Point usage | OK |
| Upper triangular ordering | OK |
| Face vertices | OK |
| Number of regions | 1 (OK) |
| Patch surfaces | each singly connected; all patches together closed |

## Table 2: Geometry and quality (checkMesh)

| Metric | Value | checkMesh limit | Status |
|---|---|---|---|
| Bounding box | (−8.059, 0, −8.059) to (16.118, 6.447, 8.059) m | — | matches `.geo` box |
| Total volume | 2512.342789 m³ | — | matches CAD fluid volume |
| Cell volume | min 1.21e−9, max 1.67 m³ | > 0 | OK |
| Face area | min 1.34e−6, max 2.94 m² | > 0 | OK |
| Boundary / cell openness | ~1e−16 | ~0 | OK |
| Max aspect ratio | 8.98 | 1000 | OK |
| ★ Non-orthogonality | max **69.3°**, avg 22.0° | 70° | OK, near limit |
| Max skewness | 1.26 | 4 | OK |
| Face pyramids | positive | — | OK |

## ★ Items for review

1. ★ **Non-orthogonality max 69.3°** passes but is close to the 70° limit. The solver may need non-orthogonal correctors.
   Location of the worst cells not yet checked (`checkMesh -writeFields '(nonOrthoAngle)'`).
2. ★ **Patch types are all `patch`** (gmshToFoam default). Before solving: `wing` → `wall`, `symmetry` → `symmetry`.
3. ★ **Tetrahedra only, no boundary layers.** Suited to inviscid (Euler) flow; a viscous (RANS) run needs prism layers on `wing`.
4. ★ **Edge refinement covers more than LE/TE.** `wingEdges()` selects 6 curves: LE line, TE line, and the root and tip
   airfoil outlines. The 3.2 mm `lc_edge` therefore also applies around the root (on the symmetry plane) and tip outlines,
   which the `.geo` comment ("fine at LE/TE") does not mention. Intended?

## Reproduce

gmsh 4.15.2; OpenFOAM v2512 (ESI/OpenCFD), Docker image `opencfd/openfoam-default:2512`.

```bash
cd week4
gmsh onera_m6_wing.geo -3 -o onera_m6_wing.msh -log gmsh.log   # ~5 s
bash run_foam.sh                                                 # gmshToFoam + checkMesh, needs Docker running
```
