#!/usr/bin/env bash
#===============================================================================
#  Viscous stencil of chimera faces at the edges of a block face (fast tier)
#
#  Two blocks overlap in x, with chimera faces (102) whose ghost cells coincide
#  with cells of the other block (make_case.py). One evaluation of the residual
#  of a linear field (one step of a steady run, forward Euler), read back from
#  the diagnostic output:
#   - plane: 2-D, u = G y; plane3d: 3-D, u = G (y + z). The cells at the ends
#     of the chimera faces must have the y (and z) momentum residual of the
#     other cells, zero, within 1e-9 of mu G / H per unit volume: with the
#     ghost cell copied beyond the edge of the block face the chimera face would
#     take 3/4 of the tangential derivative there, a residual of 1/4 of it;
#   - wedge: a 2-D axisymmetric wedge of 1 degree, a radial velocity
#     proportional to the radius: every cell of the two blocks must have the
#     residual of the same cell in one block (wedgeref), within 1e-9 of
#     mu G / H: with the ghost cell itself as the k neighbours of the ghost
#     row, the chimera faces would take half of the v/r of the divergence.
#  And a control: fanchim, plane with the y spacing of block 1 growing 4
#  times across its last cell, which turns over the first layer of ghost
#  cells extrapolated behind the chimera face, whose viscous flux takes
#  their metric: the run must stop at set-up, with a non-zero exit code,
#  on the ghost cell behind face 2.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/* OUTPUT/diagnostic-*.tec
for f in phase.txt thermo.dat transport.dat; do
  ln -sf "../../../../common/Air-viscous/$f" "INPUT/$f"
done

cat > input.ini <<'EOF'
[MOSE-Parameters]
iter-threshold = 1
newrun = true

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
EOF

status=0
for c in plane plane3d wedge wedgeref; do
  python3 -B make_case.py case "$c" || exit 1
  rm -f OUTPUT/diagnostic.tec
  OMP_NUM_THREADS=1 "$LOCAL" > "OUTPUT/$c.log" 2>&1
  rc=$?
  if [[ $rc -ne 0 || ! -f OUTPUT/diagnostic.tec ]]; then
    echo "ChimeraFace: FAIL — $c: exit code $rc, no diagnostic output (see OUTPUT/$c.log)"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$c.tec"
done

echo "-- residual of a linear field next to the edges of chimera faces --"
python3 -B make_case.py check || status=1
echo "-- a turned-over first layer of ghost cells behind a chimera face: stop at set-up --"
python3 -B make_case.py case fanchim || exit 1
rm -f OUTPUT/diagnostic.tec
OMP_NUM_THREADS=1 "$LOCAL" > OUTPUT/fanchim.log 2>&1
rc=$?
if [[ $rc -ne 0 && ! -f OUTPUT/diagnostic.tec ]] && python3 -B make_case.py stopped OUTPUT/fanchim.log; then
  :
else
  echo "  [FAIL] fanchim: exit code $rc, no stop at set-up on the ghost cell behind face 2 (see OUTPUT/fanchim.log)"
  status=1
fi
if [[ $status -eq 0 ]]; then
  echo "ChimeraFace: all checks passed"
else
  echo "ChimeraFace: FAILED"
fi
exit $status
