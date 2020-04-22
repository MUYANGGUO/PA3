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


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
    //get information of current processor
    int curRank, curPos[2];
    MPI_Comm_rank(comm, &curRank);
    MPI_Cart_coords(comm, curRank, 2, curPos);
    
    int totalProcessor, q
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
    MPI_Barrier(comm_row);
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_row);
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
    //functions: MPI_Comm_rank, MPI_Cart_coords, MPI_Cart_rank

    //get current&root rank and current&root position
    int curRank, rootRank, diagRank;
    int curPos[2], rootPos[2], diagPos[2];
    rootPos[0] = 0, rootPos[1] = 0;
    MPI_Comm_rank(comm, &curRank);
    MPI_Cart_coords(comm, curRank, 2, curPos);
    MPI_Cart_rank(comm, rootPos, &rootRank);

    diagPos[0] = curPos[0], diagPos[1] = curPos[0];
    MPI_Cart_rank(comm, DiagCoord, &DiagRank);

    //splite mpi_comm into rows and cols
    MPI_Comm comm_row, comm_col;
    MPI_Comm_split(comm, curPos[0], curPos[1], &comm_row);
    MPI_Comm_split(comm, curPos[1], curPos[0], &comm_col);

    //get R and D
    int num_row = block_decompose(n, comm_row);
    int num_col = block_decompose(n, comm_col);

    double* R = new double[num_row * num_col]();
    double* D = new double[num_col]();
    
    getRD_jacobi(R, D, num_row, num_col, loacl_A,
                curRank, diagRank, curPos, diagPos, comm_row, comm_col, comm);

    //update x until it converges or reaches to max iteration
    double* local_Rx = new double[num_col];
    double* local_Ax = new double[num_col]; 
    double l2_norm_square, sub_l2_norm_square;

    for (int i = 0; i < max_iter; i ++){
        //calculate R*x and A*x
        distributed_matrix_vector_mult(n, R, local_x, local_Rx, comm);
        distributed_matrix_vector_mult(n, local_A, local_x, local_Ax, comm);

        //update x in the first column
        if (curPos[1] == 0){
            for (int i = 0; i < num_col; i++){
                local_x[i] = (local_b[i] - local_Rx[i])/D[i];
            }
        }

        //detect termination using l2 norm
        l2_norm_square = 0;
        sub_l2_norm_square = 0;
       
        if (curPos[1] == 0){
            //get sub l2_norm of this row
            sub_l2_norm_square +=  pow(local_Ax[i] - local_b[i], 2);
        }
        MPI_Reduce(&sub_l2_norm_square, &l2_norm_square, 1, MPI_DOUBLE, MPI_SUM, 0, col_comm);

        MPI_Bcast(&l2_norm_square, 1, MPI_DOUBLE, rootRank, comm);
        if (sqrt(l2_norm_square) <= l2_termination) break;
    }

    
}

void getRD_jacobi(double* R, double* D, int num_row, int num_col, 
                double* local_A, int curRank, int diagRank, 
                int curPos[2], int diagPos[2],
                MPI_Comm comm_row, MPI_Comm comm_col, MPI_Comm comm) 
{
    //generate D and R matrix, send it to the first column
    if (curRank != diagRank){
        //get R matrix
        for (int i = 0; i < num_row * num_col; i++){
            R[i] = local_A[i];
        }
    }else{
        //get R matrix and D matrix
        for (int i = 0; i < num_row * num_col; i++){
            if (i % num_col == i / num_col){
                D[i%num_col] = local_A[i];
            }else{
                R[i] = local_A[i];
            }
        }

        //send D to the first processor in this row
        int firstProcessorRank, firstProcessorPos[2] = { curPos[0], 0};
        MPI_Cart_rank(comm, firstProcessorPos, &firstProcessorRank);
        MPI_Send(D, num_col, MPI_DOUBLE, firstProcessorRank, 1, comm);
    }

    //receive D matrix if the currant processor is at the first column
    if (curPos[1] == 0){
        MPI_Recv(D, num_col, MPI_DOUBLE, diagRank, 1, comm, MPI_STATUS_IGNORE);  
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
