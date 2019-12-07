/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>

#include "grid_base_ref_replay.h"

int main(){
    const char filename[] = "tests/QS/regtest-ot-1/grid_collocate_sample.task";
    const double max_diff = grid_collocate_replay(filename);

    if (max_diff > 1e-16) {
        printf("Max diff too high, test failed.\n");
        return 1;
    } else {
        printf("Max diff looks good, test passed.\n");
        return 0;
    }
}

//EOF
