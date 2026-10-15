"""Inputs of the WedgeAxisStep test.

A diffusion-dominated wedge, small enough that convection hardly enters the step: 1 degree wedge about the x axis, x and r
in [0, LX] x [0, LR], NX x NR cells, air of
test/common/Air-viscous at 300 K and density 1.1858686e-4 kg/m3, so nu = mu/rho = 0.1 m2/s, at rest but for a smooth radial
velocity, a swirl w = 1 m/s r/LR and a 1e-3 m/s noise; faces 300 on the i faces and on the axis, 400 on the outer j face, 200
on the k faces. With prandtl = 10 the diffusion of momentum limits the step.

With "flip" the axis is the j = NR face (r decreases with j) and the k planes are in the opposite order, so that the mesh stays
right-handed: the wedge angle of MOSE (delthe) is then negative.

Usage (from the test directory):
  python3 -B make_case.py case NX NR LX LR [flip]   INPUT/ic.tec and INPUT/bc.txt
"""
import math
import random
import sys

RGAS = 8314.46261815324 / 28.970418
T0, RHO = 300.0, 1.1858686e-4
DELTA = math.radians(1.0)


def main():
    nx, nr, lx, lr = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    flip = "flip" in sys.argv[6:]
    p0 = RHO * RGAS * T0
    dx, dr = lx / nx, lr / nr
    rn = [1.0e-6 * dr] + [j * dr for j in range(1, nr + 1)]
    if flip:
        rn = rn[::-1]                       # j = 0 the outer face, j = NR the axis
    zsign = (1.0, -1.0) if flip else (-1.0, 1.0)
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    random.seed(12345)
    u, v, w = [], [], []
    for j in range(nr):
        rc = 0.5 * (rn[j] + rn[j + 1])
        for i in range(nx):
            xc = (i + 0.5) * dx
            v.append(1.0 * (rc / lr) * math.exp(-((xc - 0.5 * lx) / (0.2 * lx)) ** 2) + 1.0e-3 * (random.random() - 0.5))
            u.append(1.0e-3 * (random.random() - 0.5))
            w.append(1.0 * rc / lr)
    with open("INPUT/ic.tec", "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        out.write(f" ZONE  T = B1-IG, I={nx + 1}, J={nr + 1}, K=2, DATAPACKING=BLOCK, "
                  "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for comp in range(3):
            for k in range(2):
                for j in range(nr + 1):
                    for i in range(nx + 1):
                        out.write(f" {(i * dx, rn[j] * c, zsign[k] * rn[j] * s)[comp]:.16E}\n")
        for col in ([RHO] * (nx * nr), u, v, w, [p0] * (nx * nr)):
            out.write("".join(f" {x:.16E}\n" for x in col))

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    axis_j, axis_face, outer_j, outer_face = (nr, 4, 1, 3) if flip else (1, 3, nr, 4)
    rows = [rec(1, j, 1, 300) for j in range(1, nr + 1)] + [rec(nx, j, 2, 300) for j in range(1, nr + 1)]
    rows += [rec(i, axis_j, axis_face, 300) for i in range(1, nx + 1)]
    rows += [rec(i, outer_j, outer_face, 400) for i in range(1, nx + 1)]
    for face in (5, 6):
        rows += [rec(i, j, face, 200) for j in range(1, nr + 1) for i in range(1, nx + 1)]
    with open("INPUT/bc.txt", "w") as out:
        out.write("".join(rows))


if __name__ == "__main__":
    if sys.argv[1] == "case":
        main()
    else:
        sys.exit(f"unknown request {sys.argv[1]}")
