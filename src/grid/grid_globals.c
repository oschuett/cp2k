/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stddef.h>
#include <stdbool.h>
#include <omp.h>

#include "grid_globals.h"

static Counters* per_thread_counters = NULL;
static bool grid_globals_initialized = false;

// *****************************************************************************
static void add_counters(const Counters increment, Counters* accumulator) {
    accumulator->collocate_ortho_cpu += increment.collocate_ortho_cpu;
    accumulator->collocate_general_cpu += increment.collocate_general_cpu;
}

// *****************************************************************************
void grid_stats_add(Counters increment) {
    assert(grid_globals_initialized);
    add_counters(increment, &per_thread_counters[omp_get_thread_num()]);
}

// *****************************************************************************
void grid_globals_init(int thread_num) {
    assert(omp_get_thread_num() == thread_num &&
           "Fortran and C runtimes disagree on OpenMP thread numbers.");

#pragma omp master
{
    assert(!grid_globals_initialized && "Grid globals were already initialized.");
    per_thread_counters = malloc(sizeof(Counters) * omp_get_num_threads());
    //TODO make this a list of pointers and allocate them in each thread individually.
    for (int i=0; i < omp_get_num_threads(); i++) {
        per_thread_counters[i] = (Counters){};
    }
    grid_globals_initialized = true;
}

#pragma omp barrier

}

// *****************************************************************************
void grid_globals_finalize(void (*mpi_sum_func)(long*), void (*print_func)(char*)){

#pragma omp barrier
#pragma omp master
{
    print_func("\n");
    print_func(" -------------------------------------------------------------------------------\n");
    print_func(" -                                                                             -\n");
    print_func(" -                                GRID STATISTICS                              -\n");
    print_func(" -                                                                             -\n");
    print_func(" -------------------------------------------------------------------------------\n");
    print_func(" COUNTER                                                                   VALUE\n");

    Counters total_counters = {};
    for (int i=0; i < omp_get_num_threads(); i++) {
        add_counters(per_thread_counters[i], &total_counters);
    }

    char buffer[100];
    mpi_sum_func(&total_counters.collocate_ortho_cpu);
    snprintf(buffer, sizeof(buffer), " %-58s %20li\n", "collocate_ortho_cpu", total_counters.collocate_ortho_cpu);
    print_func(buffer);

    mpi_sum_func(&total_counters.collocate_general_cpu);
    snprintf(buffer, sizeof(buffer), " %-58s %20li\n", "collocate_general_cpu", total_counters.collocate_general_cpu);
    print_func(buffer);

    print_func(" -------------------------------------------------------------------------------\n");

    free(per_thread_counters);
    per_thread_counters = NULL;
    grid_globals_initialized = false;
}
}

//EOF
