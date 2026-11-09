#!/usr/bin/env bash
# gmsh .msh -> OpenFOAM polyMesh -> checkMesh, inside ESI OpenFOAM v2512 (Docker).
#
# Run from anywhere:  bash week4/run_foam.sh
# Needs: Docker running, week4/onera_m6_wing.msh (from: gmsh onera_m6_wing.geo -3 -o onera_m6_wing.msh)
#
# week4/ is mounted at /work: the container sees only this folder, and writes case/constant/polyMesh here.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
IMAGE=opencfd/openfoam-default:2512

# Absolute /work paths: the image's entrypoint (/openfoam/run) cd's to /root, so -w and relative paths don't hold.
foam() { docker run --rm -v "$HERE:/work" "$IMAGE" "$@"; }

foam gmshToFoam /work/onera_m6_wing.msh -case /work/case 2>&1 | tee "$HERE/case/log.gmshToFoam"
foam checkMesh -case /work/case 2>&1 | tee "$HERE/case/log.checkMesh"
