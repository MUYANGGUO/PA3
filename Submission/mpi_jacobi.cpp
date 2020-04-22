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
    //get information of current processor
    int curRank, curPos[2];
    MPI_Comm_rank(comm, &curRank);
    MPI_Cart_coords(comm, curRank, 2, curPos);

    int totalProcessor, q;
    MPI_Comm_size(comm, &totalProcessor);
    q = (int) sqrt(totalProcessor);
    int cur_totalRows = block_decompose(n, q, curPos[0]);
    int cur_totalCols = block_decompose(n, q, curPos[1]);


    //First  step - distribute data to the first cloumn
    MPI_Comm comm_col;
    MPI_Comm_split(comm, curPos[1], curPos[0], &comm_col);

    int *row_sc = new int [q], *row_disp = new int[q];
    for (int i = 0; i < q; i++){
        row_sc[i] = n * block_decompose(n, q, i);
        row_disp[i] = (i == 0) ? 0 : row_sc[i-1] + row_disp[i-1];
    }

    double *rec_buffer_EntireRow = new double[cur_totalRows * n];
    if (curPos[1] == 0){

        MPI_Scatterv(input_matrix, row_sc, row_disp, MPI_DOUBLE,
                     rec_buffer_EntireRow, cur_totalRows * n, MPI_DOUBLE, 0, comm_col);
        MPI_Barrier(comm_col);
    }

    //Second step - distribute data to the entire row
    MPI_Comm comm_row;
    MPI_Comm_split(comm, curPos[0], curPos[1], &comm_row);

    int *col_sc = new int[q], *col_disp = new int[q];
    for (int i = 0; i < q; i++){
        col_sc[i] = block_decompose(n, q, i);
        col_disp[i] = (i == 0) ? 0 : col_sc[i-1] + col_disp[i-1];
    }

    double *rec_buffer_curP = new double[cur_totalCols * cur_totalRows];
    for (int i = 0; i < cur_totalRows; i++){
        MPI_Scatterv(rec_buffer_EntireRow + i * n, col_sc, col_disp, MPI_DOUBLE,
                    rec_buffer_curP + i * cur_totalCols, cur_totalCols, MPI_DOUBLE, 0, comm_row);

    }
    *local_matrix = rec_buffer_curP;

    MPI_Barrier(comm_row);
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_row);
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    int q = sqrt(size);
    //coords: i,j
    int coord_i = rank / q;
    int coord_j = rank % q;
    int local_m = block_decompose(n, q, coord_i);
    int local_n = block_decompose(n, q, coord_j);

    // send local data to diagonal processors from first column
    // processor(0,0)
    if (coordinates[0] == 0 && coordinates[1] == 0)
    {
        for (int i = 0; i < num_rows; i++) {
            row_vector[i] = col_vector[i];
        }
    }
    // send to diagonals
    else if (coord_j == 0)
    {
        int send_coords[] = {coord_i, coord_i};
        int send_rank;
        MPI_Cart_rank(comm, send_coords, &send_rank);
        MPI_Send(&col_vector[0], local_m, MPI_DOUBLE, send_rank, 1, comm);
    }
    // disgonals waiting to receive
    else if (coord_i == coord_j)
    {
        int original_coordinates[] = {coordinates[0], 0};
        int original_rank;
        MPI_Cart_rank(comm, original_coordinates, &original_rank);
        MPI_Recv(&row_vector[0], num_rows, MPI_DOUBLE, original_rank, 1, comm, MPI_STATUS_IGNORE);
    }

    //Bcast in each column
    MPI_Comm cols;
    MPI_Comm_split(comm, coord_j, coord_i, &cols);
    MPI_Bcast(row_vector, local_n, MPI_DOUBLE, coord_j, cols);

    MPI_Comm_free(&cols);

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
    int local_rows, local_cols;
    local_rows = block_decompose(n,rows);
    local_cols = block_decompose(n,cols);

    //transpose the local_x from (0,1) to every col
    double *x_t;
    x_t = (double*)malloc(sizeof(double)*local_rows);
    transpose_bcast_vector(n,local_x, x_t,comm);

    // cout<<"This is "<< coord_i<<" "<< coord_j<< " , local_x "<<  " is "<< *(local_x+0)<<endl;
    // if (coord_i == 0 && coord_j == 1){
    //     for (int i = 0; i < 4; i++){
    //         cout<<"This is "<< coord_i<<" "<< coord_j<< " , local_A "<< i << " is "<< *(local_A+i)<<endl;
    //     }
    // }

    //after transposing, we then could calculate the y = A*x, denote as local_nv;
    double *Ax;
    Ax = (double*)malloc(sizeof(double)*local_rows);
    for (int i = 0; i<local_rows;i++){
        Ax[i] = 0;
        for (int j = 0; j<local_cols;j++){
            Ax[i]+=local_A[ local_cols * i + j] * x_t[j];
        }
    }

    //calculate local_y, by reducing
    MPI_Reduce(Ax, local_y, local_rows, MPI_DOUBLE,MPI_SUM,0,rows);
}

void getRD_jacobi(double* R, double* D, int num_row, int num_col,
                double* local_A, int curRank, int diagRank,
                int curPos[2], int diagPos[2],
                MPI_Comm comm_row, MPI_Comm comm_col, MPI_Comm comm)
{
    //generate D and R matrix, send it to the first column
    if (curPos[0] != curPos[1]){
        //get R matrix
        for (int i = 0; i < num_row; i++){
            for (int j = 0; j < num_col; j++){
                R[i*num_col + j] = local_A[i * num_col + j];
            }
        }
    }else{
        //get R matrix and D matrix
        for (int i = 0; i < num_row; i++){
            for (int j = 0; j < num_col; j++){
                if (i == j){
                    D[i] = local_A[i * num_col + j];
                    R[i * num_col + j] = 0;
                }else{
                    R[i * num_col + j] = local_A[i * num_col + j];
                }
            }
        }

        //send D to the first processor in this row
        int firstProcessorRank, firstProcessorPos[2] = { curPos[0], 0};
        MPI_Cart_rank(comm, firstProcessorPos, &firstProcessorRank);
        if (curPos[0] != 0){
            MPI_Send(D, num_col, MPI_DOUBLE, firstProcessorRank, 1, comm);
        }
    }

    //receive D matrix if the currant processor is at the first column
    if (curPos[1] == 0 && curPos[0] != 0){
        MPI_Recv(D, num_col, MPI_DOUBLE, diagRank, 1, comm, MPI_STATUS_IGNORE);
    }
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    //functions: MPI_Comm_rank, MPI_Cart_coords, MPI_Cart_rank

    //get current&root rank and current&root position
    int curRank, rootRank, diagRank;
    int curPos[2], rootPos[2], diagPos[2];
    rootPos[0] = 0, rootPos[1] = 0;
    MPI_Comm_rank(comm, &curRank);
    MPI_Cart_coords(comm, curRank, 2, curPos);
    MPI_Cart_rank(comm, rootPos, &rootRank);

    diagPos[0] = curPos[0], diagPos[1] = curPos[0];
    MPI_Cart_rank(comm, diagPos, &diagRank);

    //splite mpi_comm into rows and cols
    MPI_Comm comm_row, comm_col;
    MPI_Comm_split(comm, curPos[0], curPos[1], &comm_row);
    MPI_Comm_split(comm, curPos[1], curPos[0], &comm_col);

    //get R and D
    int num_row = block_decompose(n, comm_row);
    int num_col = block_decompose(n, comm_col);

    double* R = new double[num_row * num_col];
    double* D = new double[num_col];

    getRD_jacobi(R, D, num_row, num_col, local_A,
                curRank, diagRank, curPos, diagPos, comm_row, comm_col, comm);

    //update x until it converges or reaches to max iteration
    double* local_Rx = new double[num_row];
    double* local_Ax = new double[num_row];
    double l2_norm_square, sub_l2_norm_square;

    // initialization
    for (int i = 0; i < num_row; i++){
        local_x[i] = 1;
    }

    for (int it = 0; it < max_iter; it++){//
        //calculate R*x and A*x
        distributed_matrix_vector_mult(n, R, local_x, local_Rx, comm);
        distributed_matrix_vector_mult(n, local_A, local_x, local_Ax, comm);
        // if (curPos[1] == 0){
        //     for (int i = 0; i < num_col; i++){
        //         cout<<"The local "<<i <<" num is "<<local_Ax[i]<<endl;
        //     }
        // }
        //detect termination using l2 norm
        l2_norm_square = 0;
        sub_l2_norm_square = 0;

        if (curPos[1] == 0){
            //get sub l2_norm of this row
            for (int i = 0; i < num_row; i++){
                sub_l2_norm_square +=  (local_Ax[i] - local_b[i])*(local_Ax[i] - local_b[i]);
            }
        }
        MPI_Reduce(&sub_l2_norm_square, &l2_norm_square, 1, MPI_DOUBLE, MPI_SUM, 0, comm_col);

        MPI_Bcast(&l2_norm_square, 1, MPI_DOUBLE, rootRank, comm);
        if (sqrt(l2_norm_square) <= l2_termination) break;
        else if (curPos[1] == 0){
            for (int i = 0; i < num_col; i++){
                local_x[i] = (local_b[i] - local_Rx[i])/D[i];
            }
        }

    }
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
