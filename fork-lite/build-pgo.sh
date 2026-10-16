#!/bin/bash
# fork-lite: profile-guided build of ngspice for a given machine.
#
# Builds twice: once with -fprofile-generate, runs the training netlists,
# then rebuilds with -fprofile-use. The result is the same code with the
# same arithmetic (no fast-math), only laid out and inlined according to
# the measured branch and call frequencies. On a Ryzen 7 3700U it was
# 7-27 % faster than the plain -O3 build on the circuits it was trained on.
#
# OpenMP is left out on purpose: on circuits below a few hundred
# transistors the parallel device loops lose to the serial path (see
# commit "devices: do not go parallel below omp_min_inst"), and the
# serial build also skips the per-instance staging of matrix entries.
#
# Usage:  fork-lite/build-pgo.sh <prefix> <training netlist>...
#   e.g.  fork-lite/build-pgo.sh /headless/ngspice47pgo tb1.spice tb2.spice
# Run the training netlists in batch mode from their own directory; they
# should resemble the circuits you simulate most (PDK, device models,
# behavioural sources).
set -e
PREFIX=${1:?prefix}; shift
[ $# -ge 1 ] || { echo "need at least one training netlist"; exit 1; }
SRC=$(cd "$(dirname "$0")/.." && pwd)
ARCH=${ARCH:--march=native}
PGODIR=${PGODIR:-/tmp/pgo-data}
CONF="--with-x=no --with-readline=yes --enable-klu --enable-osdi --enable-xspice --enable-cider --enable-predictor --disable-debug --disable-openmp"

[ -x "$SRC/configure" ] || (cd "$SRC" && ./autogen.sh)
rm -rf "$PGODIR" "$SRC/build-pgo"; mkdir -p "$PGODIR" "$SRC/build-pgo"

echo "== 1/3 instrumented build"
cd "$SRC/build-pgo"
GEN="-fprofile-generate -fprofile-dir=$PGODIR -fprofile-update=atomic"
CFLAGS="-O3 $ARCH -pipe $GEN" LDFLAGS="$GEN" ../configure --prefix="$PREFIX-gen" $CONF > configure-gen.log
make -j"$(nproc)" > make-gen.log && make install > /dev/null
# configure's own conftest programs leave profiles behind that break the
# second configure (-Werror=coverage-mismatch)
find "$PGODIR" -name '*conftest*' -delete

echo "== 2/3 training"
for n in "$@"; do
  ( cd "$(dirname "$n")" && OMP_NUM_THREADS=1 "$PREFIX-gen/bin/ngspice" -b "$(basename "$n")" > /dev/null 2>&1 ) && echo "  $n"
done
find "$PGODIR" -name '*conftest*' -delete
echo "  $(find "$PGODIR" -name '*.gcda' | wc -l) profile files"

echo "== 3/3 optimised build"
make distclean > /dev/null
USE="-fprofile-use -fprofile-dir=$PGODIR -fprofile-correction -Wno-missing-profile -Wno-error=coverage-mismatch"
CFLAGS="-O3 $ARCH -pipe $USE" LDFLAGS="$USE" ../configure --prefix="$PREFIX" $CONF > configure-use.log
make -j"$(nproc)" > make-use.log && make install > /dev/null
rm -rf "$PREFIX-gen"
echo "done: $PREFIX/bin/ngspice"
