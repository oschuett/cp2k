/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>

#include "grid_base_ref_replay.h"

int main(int argc, char *argv[]){
    if (argc != 2) {
        printf("Usage: grid_base_ref_miniapp.x <task-file>\n");
        return 1;
    }
    const int cycles = 1000;  // For better statistics the task is collocated many times.
    const double max_diff = grid_collocate_replay(argv[1], cycles);
    assert(max_diff < 1e-16);
    return 0;
}

//EOF
