/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"
#include "utils.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <cstdlib>
#include <vector>


// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
	for (int i=0; i<n; i++){
		y[i] = 0;
		for (int j=0; j<n; j++)
			y[i] += (A[i*n+j] * x[j]);
	}
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
	for (int i=0; i<n; i++){
		y[i] = 0;
		for (int j=0; j<m; j++)
			y[i] += (A[i*m+j] * x[j]);
	}
}



//function to update x
void update_x(const double* invD, const double* b, const double* R, double* x, const int n)
{
	double* y = new double[n]();
	matrix_vector_mult(n, R, x, y);
	vec_subtract(y, b, n);
	for(int i=0; i<n; i++)
		x[i] = invD[i]*y[i];
	delete[] y;

}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
	int dim = n*n;
	int iter = 0;

	double* invD = new double[n]();
	double* R = new double[dim]();
	double* y = new double[n]();

	// initialize x = [0, 0..., 0]
	for (int i=0; i<n; i++)
		x[i] = 0;

	cal_DR(A, n, invD, R);

	matrix_vector_mult(n, A, x, y);
	vec_subtract(y, b, n);
	
	while (L2_norm(n, y) > l2_termination && iter < max_iter)
	{
		iter++;		
		update_x(invD, b, R, x, n);
		matrix_vector_mult(n, A, x, y);
		vec_subtract(y, b, n);
	}

	delete[] invD;
	delete[] R;
	delete[] y;
}