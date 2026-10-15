#!/usr/bin/env bash
#===============================================================================
#  Patch test of the viscous fluxes next to walls of prescribed temperature
#  (fast tier): the cells with a face of type 302 on the meshes of
#  ViscousLinearPatch, with T linear, so that T_w varies along the walls, on
#  skewed cells; a velocity that is linear and zero on a no-slip wall,
#  the wall cells with nodes moved at random (momw), and the same on a wall
#  that does not turn with a rotating frame, where the velocity in the frame
#  is -Omega x r, varying along the wall (momr; also with the face metrics
#  inverse-mean and mean on parallelogram cells, at a wall 302 and at an
#  adiabatic wall 301, where the energy of every cell must be the work of
#  the uniform stress, the work of the wall stress included); and a
#  control that
#  every face metric passes: parallelogram cells, T_w uniform along the
#  walls (paraT). T linear at uniform pressure, the density not uniform, next
#  to walls 301, 302 and 400 faces: their ghost cells extrapolate T itself
#  (enep); a cold 302 wall and a steep pressure at 400 faces, where the floor
#  and the cut of those ghost cells act and are counted (hot). Same case
#  generator and checks: ../viscous-linear-patch
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATCH_TEST_NAME=ViscousWallPatch
export PATCH_TEST_DIR="$HERE"
export PATCH_TEST_MAKE_CASE="$HERE/../viscous-linear-patch/make_case.py"
exec bash "$HERE/../viscous-linear-patch/run.sh"
