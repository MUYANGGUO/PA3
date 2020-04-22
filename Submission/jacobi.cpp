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

    //step 1. initialize x vector, and get Diagonal matrix's inverse, R matrix
    double D[n*n];
    double R[n*n];
    for (int i = 0; i<n;i++){
        x[i] = 0;
        for (int j = 0; j<n; j++){
            if(i==j){
                //capture diagonal when col = row, and same index, R is 0
                D[i*n+j] = 1./A[i*n+j]; // inverse for later use
                R[i*n+j] = 0;
            }
            else{
                D[i*n+j] = 0;
                R[i*n+j] = A[i*n+j];
            }
        }
    }

    //step 2. while loop
    int iter = 0;
    double norm = l2_termination+1; //ensure the enter of the loop
    while (norm>l2_termination && iter<max_iter){
    
        //1.calculate || Ax - b || norm:
        //  initial a y vector by Ax
        double y[n];
        matrix_vector_mult(n,A,x,y);
        //  calculate Ax-b as norm_v
        double norm_v[n];
        for (int k =0;k<n; k++){
            norm_v[k] = y[k] - b[k];
        }
        //  calculate ||norm_v||
        for (int i=0; i<n;i++){
            norm +=norm_v[i]*norm_v[i];
        }
        norm = sqrt(norm);

        //2. update x
        //  calculate Rx
        double Rx[n];
        matrix_vector_mult(n,R,x,Rx);
        //  calculate b-Rx
        double b_Rx[n];
        for (int i = 0; i<n;i++){
            b_Rx[i] = b[i]-Rx[i];
        }
        // calculate D-1(b-Rx) to update x, here D-1 noted as D
        // update x
        matrix_vector_mult(n,D,b_Rx,x);
        iter++;
        


    }

    


    
}