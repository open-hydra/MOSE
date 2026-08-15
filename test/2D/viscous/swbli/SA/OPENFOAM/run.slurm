#!/bin/bash
#SBATCH --partition=premium
#SBATCH --job-name=swbli-saN
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ---------------------------------------------------------------
# SWBLI - rhoCentralFoam (SA, LTS) on monolith (OpenFOAM-10, .org)
#
# OpenFOAM environment comes from the login shell (~/.bashrc).
# The solver is linked against OpenMPI 4.1.4, which has no SLURM
# PMI support -> launch with its own mpirun, NOT srun.
# NB: ntasks must match numberOfSubdomains in system/decomposeParDict.
# ---------------------------------------------------------------

MPIRUN=/opt/share/comps/gcc-12.2.1/ompi-4.1.4_nccl/bin/mpirun

cd "$SLURM_SUBMIT_DIR"

# Decompose here, on the compute node, rather than on the login node: a
# login-node decomposePar dies with the ssh session and leaves a partial
# (silently broken) decomposition behind. Redo it whenever the processor
# count does not match numberOfSubdomains.
NEEDED=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
HAVE=$(ls -d processor* 2>/dev/null | wc -l)
if [ "$HAVE" -ne "$NEEDED" ]; then
    rm -rf processor*
    decomposePar -force > log.decomposePar 2>&1 || exit 1
fi

"$MPIRUN" -np "$SLURM_NTASKS" rhoCentralFoam -parallel > log.rhoCentralFoam 2>&1

reconstructPar -latestTime > log.reconstructPar 2>&1
