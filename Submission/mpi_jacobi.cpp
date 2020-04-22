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
        sendcounts = (int*)malloc(sizeof(int)*q);
        displs=(int*)malloc(sizeof(int)*q);
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
        int *recvcounts;
        int *displs;
        recvcounts = (int*)malloc(sizeof(int)*q);
        displs=(int*)malloc(sizeof(int)*q);
        //map the recvcouts, displs;
        for (int i = 0; i<q;i++){
            recvcounts[i] = block_decompose(n,q,i);
            displs[i] = (i==0)? 0:(displs[i-1]+recvcounts[i-1]);
        }

        MPI_Gatherv(local_vector,local_size,MPI_DOUBLE,output_vector,recvcounts,displs,MPI_DOUBLE,0,col);
        free(recvcounts);
        free(displs);


        
    }
    MPI_Comm_free(&col);
    return;    

}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    //get col, row from cart
    MPI_Comm row, col;
    int remain_dims_row[2] = {false, true};
    int remain_dims_col[2] = {true, false};
    MPI_Cart_sub(comm,remain_dims_row,&row);
    MPI_Cart_sub(comm,remain_dims_col,&col);




    int q = sqrt(size);
    //coords: i,j
    int coord_i = rank / q;
    int coord_j = rank % q;
    //get local size, use block_decompose function
    int local_size;
    local_size = block_decompose(n,q,coord_i);


    //send col_vector to processor at each column's diagonal position
    //meaning from each in first column (i,0), send to (i,i) respectively.
    //And diagonal processors from (i,i) should receive it.
    if (coord_j == 0){
        int Scount = block_decompose_by_dim(n,comm,0);
        MPI_Send(col_vector,Scount, MPI_DOUBLE,coord_i,coord_i,row);

    }
    if (coord_j == coord_i){
        int Rcount = block_decompose_by_dim(n,comm,0);
        MPI_Recv(row_vector,Rcount,MPI_DOUBLE,0,coord_i,row, MPI_STATUS_IGNORE);
    }
    //sync this step, wait all finished send/recv
    MPI_Barrier(comm);
    //Start Broadcast
    int Bcount = block_decompose_by_dim(n,comm,1);
    MPI_Bcast(row_vector,Bcount,MPI_DOUBLE,coord_j,col);
    //sync again after Bcast before free the memory
    MPI_Barrier(comm);

    //free the memory
    MPI_Comm_free(&col);
    MPI_Comm_free(&row);
    return;    

}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    
    int q = sqrt(size);
    //coords: i,j
    int coord_i = rank / q;
    int coord_j = rank % q;

    //MPI Split Comm to rows and cols
    MPI_Comm rows, cols;
    MPI_Comm_split(comm, coord_i,coord_j,&rows);
    MPI_Comm_split(comm, coord_j,coord_i,&cols);

    //calculate local matrix, m x n in terms of size and allocate memory acccordingly
    int local_m, local_n;
    local_m = block_decompose(n,rows);
    local_n = block_decompose(n,cols);

    //transpose the local_x from (i,0) to every row
    double *local_mv;
    local_mv = (double*)malloc(sizeof(double)*local_m);
    transpose_bcast_vector(n,local_x,local_mv,comm);

    //after transposing, we then could calculate the y = A*x, denote as local_nv;
    double *local_nv;
    local_nv = (double*)malloc(sizeof(double)*local_n);
    for (int i = 0; i<local_n;i++){
        for (int j = 0; j<local_m;j++){
            local_nv[i]+=local_A[local_n*i+j]*local_mv[j];
        }
    }

    //calculate local_y, by reducing
    MPI_Reduce(&local_nv[0],&local_y[0],local_n,MPI_DOUBLE,MPI_SUM,0,rows);
    free(local_nv);
    free(local_mv);

    MPI_Comm_free(&rows);
    MPI_Comm_free(&cols);


    return;

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
