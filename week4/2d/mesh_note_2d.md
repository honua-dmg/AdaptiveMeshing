# ONERA D airfoil (2-D): mesh

**2-D mesh: 47,449 cells (25,525 triangles + 21,924 boundary-layer quads), 34,951 nodes.**
**Extruded mesh: 47,449 cells (25,525 prisms + 21,924 hexahedra), one cell thick. checkMesh: 1 failed check (★1, by design).**

Generated from the supplied `onera_d_airfoil_2d.geo` (`VISCOUS = 1`, quad boundary layer, first cell 1.571 µm),
meshed in 2-D (`EXTRUDE = 0`) and extruded in OpenFOAM with `extrudeMesh`, because the `.geo`'s own extrusion gives
broken cells at the wall (★2). Items marked ★ need review.

## Files

| File | What |
|---|---|
| `onera_d_airfoil_2d.msh` | 2-D gmsh mesh, format 2.2 ASCII, 3.3 MB |
| `case/constant/polyMesh/` | extruded OpenFOAM mesh |
| `gmsh.log`, `case/log.extrudeMesh`, `case/log.topoSet`, `case/log.createPatch`, `case/log.checkMesh` | full tool output |

| Patch | 2-D mesh (edges) | Extruded mesh (faces) | OpenFOAM type |
|---|---|---|---|
| `airfoil` | 424 | 424 | `wall` |
| `farfield` | 105 | 105 | `patch` |
| `front` | — | 47,449 | `empty` |
| `back` | — | 47,449 | `empty` |

## Table 1: Topology (checkMesh, extruded mesh)

| Check | Result |
|---|---|
| Boundary definition | OK |
| Cell-to-face addressing | OK |
| Point usage | OK |
| Upper triangular ordering | OK |
| Face vertices | OK |
| Number of regions | 1 (OK) |
| Patch surfaces | each singly connected; all patches together closed |

## Table 2: Geometry and quality (checkMesh, extruded mesh)

| Metric | Value | checkMesh limit | Status |
|---|---|---|---|
| Bounding box | (−40.08, −40.29, 0) to (40.50, 40.29, 0.0806) m | — | far-field circle R = 50c, span 0.1c |
| Solution directions | 2, (1 1 0) | — | 2-D case recognised |
| Edges aligned with 2-D plane | all | — | OK |
| Total volume | 410.8377 m³ | — | = polygonal far-field area × span |
| Cell volume | min 4.59e−14, max 0.283 m³ | > 0 | OK |
| Face area | min 5.70e−13, max 3.52 m² | > 0 | OK |
| Boundary openness | ~1e−15 | ~0 | OK |
| ★ Max aspect ratio | **1406** (8 cells) | 1000 | **failed** |
| Non-orthogonality | max 55.3°, avg 3.4° | 70° | OK |
| Max skewness | 0.61 | 4 | OK |
| Face pyramids | positive | — | OK |

## ★ Items for review

1. ★ **Aspect ratio 1406 in 8 cells** (checkMesh limit 1000). Expected for a y+ ~ 1 boundary layer: first cell 1.571 µm
   normal to the wall, up to ~2 mm along it. Wall-aligned stretched cells are normally acceptable; the 8 cells are
   written to the set `highAspectRatioCells` (location not yet inspected).
2. ★ **The `.geo`'s own extrusion (`EXTRUDE = 1`) gives broken cells at the wall** with gmsh 4.15.2: 1,075 hexahedra,
   6 prisms, 900 quads and 6 triangles with repeated nodes (zero volume) along the airfoil boundary layer.
   `gmshToFoam` then aborts (`Patch front starts at 80331 / Current face counter at 80175`). Meshing with
   `EXTRUDE = 0` (0 collapsed elements) or `VISCOUS = 0` (0 collapsed) is clean: only the combination fails.
   Not tested with other gmsh versions.
3. ★ **`extrudeMesh` default point merging would destroy the boundary layer.** It merges points closer than
   1e−4 × domain size (8.06 µm here), which collapsed 10,392 edges of the 1.571 µm wall cells when tested.
   Set `mergeTol 1e-12` in `system/extrudeMeshDict`; 0 edges merged.
4. ★ **Patch split by position.** `extrudeMesh` puts all side faces in one patch; `topoSet` + `createPatch` assign faces
   inside the box (−0.1, −0.1)–(0.9, 0.1) to `airfoil` and the rest to `farfield`. Counts match the `.geo` (424 / 105).

## Reproduce

gmsh 4.15.2; OpenFOAM v2512 (ESI/OpenCFD), Docker image `opencfd/openfoam-default:2512`.

```bash
cd week4/2d
GMSH=gmsh bash run_foam_2d.sh     # 2-D mesh (EXTRUDE = 0 copy of the .geo) -> extrudeMesh -> topoSet -> createPatch -> checkMesh
```
