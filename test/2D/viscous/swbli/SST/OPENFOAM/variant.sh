#!/bin/bash
# Select the SST variant of this case. The names are those of the reference
# wall solutions, ../reference/OF-SST-<wall>-<production>.xy:
#
#   <wall>        asymptotic    MOSE's omega wall condition: 80 nu_w/y_c^2 on the
#                               wall FACE, wall-adjacent cell solved (system/omegaFaceBC)
#                 wallfunction  OpenFOAM's omegaWallFunction: the same value FIXED
#                               in the wall-adjacent cell
#   <production>  compressible   kOmegaSST, P = tau_ij du_i/dx_j
#                                (= MOSE sst-production = compressible)
#                 incompressible kOmegaSSTMOSE, P = mu_t S^2
#                                (= MOSE sst-production = incompressible, the default);
#                                built here from kOmegaSSTMOSE/ if not yet compiled
#
# usage: ./variant.sh <asymptotic|wallfunction> <compressible|incompressible>
# The committed case is "asymptotic compressible".
cd "${0%/*}" || exit 1
set -e

wall=$1; prod=$2
case "$wall" in asymptotic|wallfunction) ;; *) echo "usage: $0 <asymptotic|wallfunction> <compressible|incompressible>"; exit 1;; esac
case "$prod" in compressible|incompressible) ;; *) echo "usage: $0 <asymptotic|wallfunction> <compressible|incompressible>"; exit 1;; esac

# omega on bottomWall and topWall (the inlet value line is left alone)
if [ "$wall" = wallfunction ]; then
    sed -i '/#include "$FOAM_CASE\/system\/omegaFaceBC"/{n;s|uniform 1e9;|uniform 1215.8;|}' 0/omega
    sed -i 's|#include "$FOAM_CASE/system/omegaFaceBC"|type            omegaWallFunction;|' 0/omega
else
    sed -i '/type            omegaWallFunction;/{n;s|uniform 1215.8;|uniform 1e9;|}' 0/omega
    sed -i 's|type            omegaWallFunction;|#include "$FOAM_CASE/system/omegaFaceBC"|' 0/omega
fi

# production: model and library
sed -i '/^libs .*libkOmegaSSTMOSE/d' system/controlDict
if [ "$prod" = incompressible ]; then
    [ -f "$FOAM_USER_LIBBIN/libkOmegaSSTMOSE.so" ] || wmake libso kOmegaSSTMOSE
    sed -i 's/^\( *model *\)kOmegaSST\(MOSE\)\?;/\1kOmegaSSTMOSE;/' constant/momentumTransport
    sed -i '/^application/a libs            ("libkOmegaSSTMOSE.so");' system/controlDict
else
    sed -i 's/^\( *model *\)kOmegaSST\(MOSE\)\?;/\1kOmegaSST;/' constant/momentumTransport
fi

echo "variant: $wall $prod"
grep -E '^ *(#include .*omegaFaceBC"|type +omegaWallFunction)' 0/omega | sed 's/^ */  0\/omega: /'
grep -E '^ *model' constant/momentumTransport | sed 's/^ */  momentumTransport: /'
grep -E '^libs' system/controlDict | sed 's/^/  controlDict: /' || true
