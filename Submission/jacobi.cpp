/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>
using namespace std;

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    // TODO
    for (int i = 0; i<n; i++){
        y[i] = 0;

        for (int j=0; j<n; j++){
            y[i] += (A[i*n+j]*x[j]);
        }
    }
    //   cout<< "\nHello World!"<<*y<<"check\n";
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
    for (int i = 0; i<n; i++){
        y[i] = 0;

        for (int j=0; j<m; j++){
            y[i] += (A[i*m+j]*x[j]);
        }
    }

}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    // TODO
}
