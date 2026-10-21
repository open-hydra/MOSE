#!/usr/bin/env bash
#===============================================================================
#  Restart round-trip (fast tier, I3/F6)
#
#  Verifies that stopping a run, writing a solution and restarting (newrun=false,
#  which reads OUTPUT/field.tec) reproduces the uninterrupted run:
#     A : N iterations in one go
#     B : N/2 iterations -> restart -> N/2 more
#  A and B must match to solution-file round-trip precision.  Reuses the coarse
#  Sod inputs (test/1D/Sod79); no new mesh/IC data is stored here.
#     C : probes at a restart. Two probes on cells of the tube, by index (pA on
#         cell 50, pB on cell 45: the nodes of this 1-D mesh lie on the planes
#         j = 0 and k = 0), one row per iteration: 20 iterations in one go,
#         then 10 with the file of pA removed before a restart for 10 more.
#         The restart must end normally, create pA.txt again with the rows of
#         the restarted part and append those rows to pB.txt, and the rows of
#         both files must be those of the run in one go (check_probes.py; the
#         restarted run numbers its iterations from 1 again). A
#         restart did not open the probes from a missing file on, and the run
#         ended at its first probe output with a run-time error.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air
N=50          # total iterations
HALF=25       # iterations before the restart

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt

write_ini() {  # <iter-threshold> <newrun>
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = $2

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.9
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
}

# --- A: uninterrupted N-iteration run --------------------------------------
write_ini "$N" true
rm -f OUTPUT/field*.tec
"$LOCAL" > OUTPUT/solve_full.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from full run"; exit 1; }
cp OUTPUT/field.tec field_full.tec

# --- B: N/2 then restart for the remaining N/2 -----------------------------
write_ini "$HALF" true
rm -f OUTPUT/field*.tec
"$LOCAL" > OUTPUT/solve_half.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from first half"; exit 1; }

write_ini "$HALF" false
"$LOCAL" > OUTPUT/solve_restart.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from restart"; exit 1; }
cp OUTPUT/field.tec field_restart.tec

python3 -B check_restart.py field_full.tec field_restart.tec
status=$?

# --- C: probes at a restart, the file of one of them removed ----------------
add_probes() {
  cat >> input.ini <<EOF

[MOSE-Probes]
probe1 = pA
probe2 = pB

[pA]
variables = p u
diter = 1
index-position = 1 50 1 1

[pB]
variables = p
diter = 1
index-position = 1 45 1 1
EOF
}
write_ini 20 true
add_probes
rm -f OUTPUT/pA.txt OUTPUT/pB.txt
"$LOCAL" > OUTPUT/probes_full.log 2>&1
mv -f OUTPUT/pA.txt OUTPUT/pA_full.txt 2>/dev/null
mv -f OUTPUT/pB.txt OUTPUT/pB_full.txt 2>/dev/null
write_ini 10 true
add_probes
rm -f OUTPUT/pA.txt OUTPUT/pB.txt
"$LOCAL" > OUTPUT/probes_first.log 2>&1
rows() { if [[ -f $1 ]]; then wc -l < "$1"; else echo 0; fi; }
nA=$(rows OUTPUT/pA.txt)
nB=$(rows OUTPUT/pB.txt)
rm -f OUTPUT/pA.txt
write_ini 10 false
add_probes
"$LOCAL" > OUTPUT/probes_restart.log 2>&1
rc=$?
mA=$(rows OUTPUT/pA.txt)
mB=$(rows OUTPUT/pB.txt)
python3 -B check_probes.py 10 OUTPUT/pA_full.txt OUTPUT/pA.txt 0 OUTPUT/pB_full.txt OUTPUT/pB.txt "$nB"
same=$?
if [[ $rc -eq 0 && $nA -gt 0 && $mA -eq $nA && $mB -eq $((2 * nB)) && $same -eq 0 ]]; then
  echo "  [ok  ] probes at a restart: pA.txt created again ($mA rows), pB.txt appended ($nB + $nB rows), the rows of the run in one go"
else
  echo "  [FAIL] probes at a restart: exit status $rc; pA.txt $mA rows (before the restart $nA), pB.txt $mB rows" \
       "(before $nB), rows of the run in one go: $([[ $same -eq 0 ]] && echo yes || echo no); see OUTPUT/probes_restart.log"
  status=1
fi
exit $status
