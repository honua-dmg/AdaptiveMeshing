#!/usr/bin/env bash
# 2-D ONERA D airfoil: gmsh 2-D mesh -> OpenFOAM extrudeMesh (one cell thick, front/back empty) -> checkMesh.
# ESI OpenFOAM v2512 in Docker.
#
# Run from anywhere:  GMSH=/path/to/gmsh bash week4/2d/run_foam_2d.sh     (GMSH defaults to "gmsh" on PATH)
# Needs: Docker running, gmsh >= 4.8.
#
# Why not the .geo's own extrusion (EXTRUDE = 1): combined with the boundary layer it produces collapsed cells
# in gmsh 4.15.2 (1,075 hexes with repeated nodes) and gmshToFoam aborts. So: mesh in 2-D, extrude in OpenFOAM.
#
# Outputs: onera_d_airfoil_2d.msh (2-D mesh), case/constant/polyMesh (extruded mesh), case/log.*
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
IMAGE=opencfd/openfoam-default:2512
GMSH="${GMSH:-gmsh}"

# 1. 2-D mesh. The .geo sets EXTRUDE = 1 by plain assignment, which overrides gmsh -setnumber, so mesh a copy
#    with EXTRUDE = 0. The requester's file is not modified.
sed 's/^EXTRUDE  = 1;/EXTRUDE  = 0;/' "$HERE/onera_d_airfoil_2d.geo" > "$HERE/_2d.geo"
trap 'rm -f "$HERE/_2d.geo"' EXIT
grep -q '^EXTRUDE  = 0;' "$HERE/_2d.geo" || { echo "EXTRUDE line not found in .geo; edit the sed above" >&2; exit 1; }
"$GMSH" "$HERE/_2d.geo" -2 -o "$HERE/onera_d_airfoil_2d.msh" -log "$HERE/gmsh.log" > /dev/null
"$GMSH" "$HERE/onera_d_airfoil_2d.msh" -0 -format obj -o "$HERE/case/surf.obj" > /dev/null   # surface for extrudeMesh

# 2. Extrude, name patches, check. Absolute /work paths: the image's entrypoint cd's to /root.
foam() { docker run --rm -v "$HERE:/work" "$IMAGE" "$@"; }
rm -rf "$HERE/case/constant"
foam extrudeMesh -case /work/case 2>&1 | tee "$HERE/case/log.extrudeMesh"
foam topoSet -case /work/case 2>&1 | tee "$HERE/case/log.topoSet"
foam createPatch -overwrite -case /work/case 2>&1 | tee "$HERE/case/log.createPatch"
foam checkMesh -case /work/case 2>&1 | tee "$HERE/case/log.checkMesh"
