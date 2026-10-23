#!/usr/bin/env bash
#===============================================================================
#  Patch test of the viscous fluxes next to inlets and outlets (fast tier):
#  every outer face of a 2-D mesh with nodes moved at random (three sets of
#  types) and of a 3-D pair of blocks an inlet or an outlet, one type per
#  block face: 401, 402, 403, 404, 405, 406, 407 and 420, whose ghost cells
#  are copies of the cell, and 408, whose ghost cells extrapolate the cells,
#  3 P1 - 3 P2 + P3. These faces have no viscous flux, so with a linear
#  velocity (mom) or a linear temperature at uniform density (ene) the cells
#  of a face must give minus the force of the uniform stress, or the uniform
#  heat conduction, through it, and the other cells zero, within 1e-11 of the
#  scale of the face fluxes: they read the ghost cells through the faces
#  between the cells of a face (in 2-D also through their symmetry faces in
#  k), at the positions that the face metric gives them. The density is
#  uniform because the 408 ghost cells extrapolate the density and the
#  pressure: with a density that varies they do not hold T on the linear
#  field. A laminar test, one rank, one thread: it does not exercise the cell
#  metrics of the RANS source terms or the MUSCL reconstruction at the first
#  interior face, whose inviscid flux drops out of the difference of the two
#  steps. Same case generator and checks: ../viscous-linear-patch
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATCH_TEST_NAME=ViscousOpenPatch
export PATCH_TEST_DIR="$HERE"
export PATCH_TEST_MAKE_CASE="$HERE/../viscous-linear-patch/make_case.py"
exec bash "$HERE/../viscous-linear-patch/run.sh"
