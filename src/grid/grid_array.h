/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_ARRAY_H
#define GRID_ARRAY_H

// *****************************************************************************
#define Array2dAt(array, i, j) array.data[i*array.s1 + j]

typedef struct {
   double* data;
   int s1; // stride 1
} Array2d;

// *****************************************************************************
#endif

//EOF