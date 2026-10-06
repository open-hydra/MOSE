#!/usr/bin/env bash
#===============================================================================
#  Patch test of the viscous fluxes next to walls of prescribed temperature
#  (fast tier): the cells with a face of type 302 on the meshes of
#  ViscousLinearPatch, with T linear, so that T_w varies along the walls, on
#  skewed cells; a velocity that is linear and zero on a no-slip wall,
#  the wall cells with nodes moved at random (momw); and a control that
#  every face metric passes: parallelogram cells, T_w uniform along the
#  walls (paraT). Same case generator and checks: ../viscous-linear-patch
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATCH_TEST_NAME=ViscousWallPatch
export PATCH_TEST_DIR="$HERE"
export PATCH_TEST_MAKE_CASE="$HERE/../viscous-linear-patch/make_case.py"
exec bash "$HERE/../viscous-linear-patch/run.sh"
