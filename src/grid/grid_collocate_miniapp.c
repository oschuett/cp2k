/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <omp.h>

#include "grid_globals.h"
#include "grid_collocate_replay.h"

void mpi_sum_func(long* number){
    *number += 0;  // Nothing todo without MPI, pretend argument is used anyways.
}

void print_func(char* message){
    printf(message);
}

int main(int argc, char *argv[]){
    if (argc != 2) {
        printf("Usage: grid_base_ref_miniapp.x <task-file>\n");
        return 1;
    }

    grid_globals_init(omp_get_thread_num());

    const int cycles = 1000;  // For better statistics the task is collocated many times.
    const double max_diff = grid_collocate_replay(argv[1], cycles);
    assert(max_diff < 1e-11 * cycles);

    grid_globals_finalize(&mpi_sum_func, &print_func);

    return 0;
}

//EOF
