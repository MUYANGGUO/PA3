/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

/*
 * TODO: Implement your solutions here
 */
#include<iostream>

using namespace std;


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO

    //Define the processors in the communicator in 2D grids (i,j)
    //root processor is root zero (0,0)
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    MPI_Comm col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm,remain_dims,&col);

    int q = sqrt(size);
    //coords: i,j
    int coord_i = rank / q;
    int coord_j = rank % q;
    //get local size, use block_decompose function
    int local_size;
    local_size = block_decompose(n,q,coord_i);
    //distribute from root (0,0), to coord_j = 0, that is (i,0) first column
    if (coord_j == 0){
        int *sendcounts;
        int *displs;
        sendcounts = (int*)malloc(sizeof(int)*size);
        displs=(int*)malloc(sizeof(int)*size);
        //check: assign memory to local vector
        *local_vector =  (double*)malloc(sizeof(double)*local_size);
        //map the sendcouts, displs;
        for (int i = 0; i<q;i++){
            sendcounts[i] = block_decompose(n,q,i);
            displs[i] = (i==0)? 0:(displs[i-1]+sendcounts[i-1]);
        }

        MPI_Scatterv(input_vector, sendcounts,displs,MPI_DOUBLE,*local_vector,local_size,MPI_DOUBLE,0,col);

        free(sendcounts);
        free(displs);


        
    }
    MPI_Comm_free(&col);
    return;

    
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
