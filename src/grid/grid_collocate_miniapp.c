/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "grid_collocate_replay.h"

int main(int argc, char *argv[]){
    bool cuda = false;
    char* filename = NULL;

    if (argc == 2) {
        filename = argv[1];
    } else if (argc == 3 && strcmp(argv[1], "--cuda") == 0) {
        cuda = true;
        filename = argv[2];
    }

    if (filename == NULL) {
        printf("Usage: grid_base_ref_miniapp.x [--cuda] <task-file>\n");
        return 1;
    }
    const int cycles = 20000;  // For better statistics the task is collocated many times.
    const double max_diff = grid_collocate_replay(filename, cycles, cuda);
    assert(max_diff < 1e-10 * cycles);
    return 0;
}

//EOF
