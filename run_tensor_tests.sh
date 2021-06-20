#!/bin/bash -e

TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-mp2-grad"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-gw-cubic"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-gw2x"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-ri-rpa"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-mp2-stress"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-hfx-ri"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-gw"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-ri-mp2"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-rpa-cubic-scaling"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-rpa-lr"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-mp2-lr"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-sos-mp2-lr"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-gw-ic-model"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-ri-rpa-rse"
TESTOPTS="${TESTOPTS} -restrictdir QS/regtest-ri-laplace-mp2-cubic"

set -x

make VERSION=sdbg test TESTOPTS="${TESTOPTS}"

#EOF
