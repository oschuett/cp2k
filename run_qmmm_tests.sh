#!/bin/bash -e

TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-1"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-nmr-6"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-2-erf"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-gapw"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-2-swave"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-3"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-shell-pol"
TESTOPTS="${TESTOPTS} -restrictdir QMMM/QS/regtest-4"

set -x

make VERSION=sdbg test TESTOPTS="${TESTOPTS}"

#make VERSION=pdbg test TESTOPTS="${TESTOPTS}"
#make ARCH=local_coverage VERSION=pdbg test TESTOPTS="${TESTOPTS}"

#EOF
