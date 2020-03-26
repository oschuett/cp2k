/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_GLOBALS_H
#define GRID_GLOBALS_H

// *****************************************************************************
void grid_globals_init(int thread_num);
void grid_globals_finalize(void (*mpi_sum_func)(long*),
                           void (*print_func)(char*));

// *****************************************************************************
typedef struct {
   long collocate_ortho_cpu;
   long collocate_general_cpu;
} Counters;

// *****************************************************************************
void grid_stats_add(Counters increment);

#endif  // GRID_GLOBALS_H

//EOF
