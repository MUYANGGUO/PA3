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
#include <iostream>
using namespace std;

/*
 * TODO: Implement your solutions here
 */


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm) {
    // TODO
    /*Equally distributes a vector stored on processor (0,0) onto
    *          processors (i,0) [i.e., the first column of the processor grid].
    *
            * Block distributes the input vector of size `n` from process (0,0) to
    * processes (i,0), i.e., the first column of the 2D grid communicator.*/

    //get rank and coord in cart comm
    int myRank, myRankCoord[2];
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, &myRankCoord[0]);

    //create group for first column
    int remain_dims[2] = {true, false};
    MPI_Comm col_comm;
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    //get p and q
    int mySize; //size of p = q*q
    MPI_Comm_size(comm, &mySize);
    int mySizeCart; // get value of q
    mySizeCart = (int) sqrt(mySize);

    int local_size = block_decompose(n, mySizeCart, myRankCoord[0]);
    if (myRankCoord[1] == 0) {
        (*local_vector) = new double[local_size]; //assign memory for local vector

        //get information for scatter
        int *sendCounts = new int[mySizeCart];
        int *disp = new int[mySizeCart];
        for (int i = 0; i < mySizeCart; i++) {
            sendCounts[i] = block_decompose(n, mySizeCart, i);
            disp[i] = (i == 0) ? 0 : (disp[i - 1] + sendCounts[i - 1]);
        }

        MPI_Scatterv(input_vector, sendCounts, disp, MPI_DOUBLE, *local_vector, local_size, MPI_DOUBLE,
                     0, col_comm);

        delete[] sendCounts;
        delete[] disp;
    }
    MPI_Comm_free(&col_comm);
    return;
}

// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    /*Reverts the operation of `distribute_vector()`.
    *
            * Gathers the vector `local_vector`, which is distributed among the first
    * column of processes (i,0), onto the process with rank (0,0).*/

    //get rank and coord in cart comm
    int myRank, myRankCoord[2];
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, &myRankCoord[0]);

    //create group for first column
    int remain_dims[2] = {true, false};
    MPI_Comm col_comm;
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    //get p and q
    int mySize; //size of p = q*q
    MPI_Comm_size(comm, &mySize);
    int mySizeCart; // get value of q
    mySizeCart = (int) sqrt(mySize);

    int local_size = block_decompose(n, mySizeCart, myRankCoord[0]);

    if (myRankCoord[1] == 0) {
        //get information for gather
        int *recvCounts = new int[mySizeCart];
        int *disp = new int[mySizeCart];
        for (int i = 0; i < mySizeCart; i++) {
            recvCounts[i] = block_decompose(n, mySizeCart, i);
            disp[i] = (i == 0) ? 0 : (disp[i - 1] + recvCounts[i - 1]);
        }

        MPI_Gatherv(local_vector,local_size,MPI_DOUBLE,output_vector,recvCounts,disp,MPI_DOUBLE,0,col_comm);

        delete[] recvCounts;
        delete[] disp;
    }
    MPI_Comm_free(&col_comm);
    return;
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
/*    Equally distributes a matrix stored on processor (0,0) onto the
    *          whole grid (i,j).
                     *
                     * Block distributes the input matrix of size n-by-n, stored in row-major
                     * format onto a 2d communicator grid of size q-by-q with a total of
                     * p = q*q processes.*/
    //get rank and coord in cart comm
    int myRank, myRankCoord[2];
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, &myRankCoord[0]);

    //create group for first column
    int remain_dims[2] = {true, false};
    MPI_Comm col_comm;
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    //get p and q
    int mySize; //size of p = q*q
    MPI_Comm_size(comm, &mySize);
    int mySizeCart; // get value of q
    mySizeCart = (int) sqrt(mySize);

    int local_size = n*block_decompose(n, mySizeCart, myRankCoord[0]);
    //if (myRankCoord[1] == 0)
        double *local_matrix_in_first_col = new double[local_size]; //assign memory for first scatter

    //get information for scatter
    int *sendCounts = new int[mySizeCart];
    int *disp = new int[mySizeCart];
    for (int i = 0; i < mySizeCart; i++) {
        sendCounts[i] = n * block_decompose(n, mySizeCart, i);
        disp[i] = (i == 0) ? 0 : (disp[i - 1] + sendCounts[i - 1]);
    }

    if (myRankCoord[1] == 0) {
        //cout << " scatter on rank [" << myRankCoord[0] << ", " << myRankCoord[1] << "]" << endl;
        MPI_Scatterv(input_matrix, sendCounts, disp, MPI_DOUBLE, local_matrix_in_first_col, local_size, MPI_DOUBLE,
                     0, col_comm);
        MPI_Barrier(col_comm);
    }

    //second step: scatter matrix from the first cols
    //create group for rows
    remain_dims[0] = false; remain_dims[1]= true;
    MPI_Comm row_comm;
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    int local_size_row = block_decompose(n, mySizeCart, myRankCoord[0]);
    int local_size_col = block_decompose(n, mySizeCart, myRankCoord[1]);
    local_size = local_size_row*local_size_col;
    (*local_matrix) = new double[local_size]; //assign memory for each block

    //get information for scatter
    int *sendCountsRow = new int[mySizeCart];
    int *dispRow = new int[mySizeCart];
    for (int i = 0; i < mySizeCart; i++) {
        sendCountsRow[i] = block_decompose(n, mySizeCart, i) ;
        dispRow[i] = (i == 0) ? 0 : (dispRow[i - 1] + sendCountsRow[i - 1]);
    }

    for (int i=0; i<local_size_row; i++) {
        MPI_Scatterv((local_matrix_in_first_col + i * n), sendCountsRow, dispRow, MPI_DOUBLE,
                     (*local_matrix + i * local_size_col), local_size_col, MPI_DOUBLE,
                     0, row_comm);
    }
    MPI_Barrier(row_comm);

    delete[] sendCounts;
    delete[] sendCountsRow;
    delete[] disp;
    delete[] dispRow;
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    return;

}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    /*Given a vector distributed among the first column,
    *          this function transposes it to be distributed among rows.
    * Given a vector that is distirbuted among the first column (i,0) of processors,
    * this function will "transpose" the vector, such that it is block decomposed
    * by row of processors.*/

    int myRank, myRankCoord[2];
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, &myRankCoord[0]);

    //get row_comm
    int remain_dims[2];
    remain_dims[0] = false; remain_dims[1]= true;
    MPI_Comm row_comm;
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    //get col_comm
    remain_dims[0] = true; remain_dims[1]= false;
    MPI_Comm col_comm;
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    //send col_vector to diag
    if(myRankCoord[1]==0){
        MPI_Send(col_vector,block_decompose_by_dim(n,comm,0),MPI_DOUBLE,myRankCoord[0],myRankCoord[0],row_comm);
    }
    if(myRankCoord[0]==myRankCoord[1]){
        MPI_Recv(row_vector,block_decompose_by_dim(n,comm,0),MPI_DOUBLE,0,myRankCoord[0],row_comm,MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    //broadcast vector in col
    MPI_Bcast(row_vector,block_decompose_by_dim(n,comm,1),MPI_DOUBLE,myRankCoord[1],col_comm);

    MPI_Barrier(comm);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);

    return;
}




void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    /** TODO
    *  Calculates y = A*x for a distributed n-by-n matrix A and
    * distributed, size n vectors x and y on a q-by-q processor grid.
    * The matrix A is distirbuted on the q-by-q grid, and the vectors are 
    * distirbuted on the q-by-q grid (block distributed accross the first column. )
    *
    * First transposing the input vector x via transpose_bcast_vector(), then locally multiply
    * the row decomposed vector by the local matrix. Then, the resulting local vectors are
    * summed by using MPI_Reduce along rows onto the first column, which yields result y.
     **/
    
    int myRank, nRow, nCol;
    int myCoords[2];
    MPI_Comm row_comm, col_comm;

    // get cart coords
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, myCoords);

    // split comm into row and column comms
    MPI_Comm_split(comm, myCoords[0], myCoords[1], &row_comm); 
    MPI_Comm_split(comm, myCoords[1], myCoords[0], &col_comm);
    
    // get block (local) row and column sizes
    nCol = block_decompose(n, col_comm);    
    nRow = block_decompose(n, row_comm);    

    // transpose local_x from first column (i,0) to each row
    // use row_x array to hold the distributed local_x vector
    double *row_x = new double[nRow];
    transpose_bcast_vector(n, local_x, row_x, comm);    

    // compute block product block y = block (A*x)
    double *sub_y= new double[nCol]();
    for (int i=0; i<nCol; i++)
        for (int j=0;j<nRow; j++)
            sub_y[i] += local_A[nRow*i+j] * row_x[j];
        
    // reduce (sum) sub_y on all processors to local_y on column zero
    MPI_Reduce(&sub_y[0], &local_y[0], nCol, MPI_DOUBLE, MPI_SUM, 0, row_comm);
    
    delete[] row_x;
    delete[] sub_y;

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
}



// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{

    int myRank, nRow, nCol, rootRank;
    int rootRankCoords[2] = {0, 0};
    int myCoords[2];

    MPI_Comm row_comm, col_comm;

    // create interchangable ranks & coords
    MPI_Comm_rank(comm, &myRank);
    MPI_Cart_coords(comm, myRank, 2, myCoords);
    MPI_Cart_rank(comm, rootRankCoords, &rootRank);

    // obtain rank and coords for diagonal elements
    // for further calculation of D and R
    int DiagRank;
    int DiagCoord[2] = {myCoords[0], myCoords[0]};
    MPI_Cart_rank(comm, DiagCoord, &DiagRank);

    MPI_Comm_split(comm, myCoords[0], myCoords[1], &row_comm); // split comm into row comms
    MPI_Comm_split(comm, myCoords[1], myCoords[0], &col_comm); // split comm into col comms
    
    // decompose row and col block sizes
    nRow = block_decompose(n, row_comm);
    nCol = block_decompose(n, col_comm);

    double* D = new double[nCol]();
    double* invD = new double[nCol]();
    double* R = new double[nRow*nCol](); 

    // initiate x vector
    for (int i=0; i<nCol; i++)
        local_x[i] = 0;

    // obtain R and D from A
    if (myRank == DiagRank){
        // cout << myCoords[0] << ", " << myCoords[1] << "\t";
        for(int i=0; i<nCol; i++){
            D[i]=local_A[i*nCol+i];
            for(int j=0; j<nCol; j++){
                if(i!=j)
                    R[i*nCol+j]= local_A[i*nCol+j];
            }
        }

        // send D to the first column
        int dest;
        int destCoords[2];
        destCoords[1] = 0;
        destCoords[0] = myCoords[0];
        MPI_Cart_rank(comm, destCoords, &dest);
        MPI_Send(D, nCol, MPI_DOUBLE, dest, 1, comm);
    }
    else
        for (int k=0; k<nRow*nCol; k++)
            R[k] = local_A[k];

    // MPI_Barrier(comm);

    // first column receives D and computes inverse(D)
    if (myCoords[1] == 0){
        int src;
        int srcCoords[2] = {myCoords[0], myCoords[0]};
        MPI_Cart_rank(comm, srcCoords, &src);
        MPI_Recv(D, nCol, MPI_DOUBLE, src, 1, comm, MPI_STATUS_IGNORE);  
        for (int i=0; i<nCol; i++){
            invD[i] = 1/D[i];
        }
    }

    double  l2_norm, temp;
    double* local_y = new double[nCol];
    double* res = new double[nCol];

    // use jacobi iteration to compute x 
    // iteration terminates when l2_norm < L2_termination
    for (int iter=0; iter<max_iter; iter++)
    {
        distributed_matrix_vector_mult(n, R, local_x, local_y, comm);
        if (myCoords[1] == 0)
        {
            for(int i=0; i<nCol; i++){
                res[i] = local_b[i] - local_y[i];
                local_x[i]=invD[i]*res[i];
            }
        }

        // distriuted matrix * vector multiplication, store results at local_y
        // compute residual and sum up
        // reduce to root [0,0] for the calculation of L2 norm
        distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);
        if (myCoords[1] == 0)
        {
            temp = 0;

            for(int i=0; i<nCol; i++){
                res[i]= local_y[i] - local_b[i];
                temp += res[i]*res[i];
            }

            MPI_Reduce(&temp, &l2_norm, 1, MPI_DOUBLE, MPI_SUM, 0, col_comm);
        }

        if (myRank == rootRank)
            l2_norm = sqrt(l2_norm);
        
        MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, rootRank, comm);

        // break if L2 Norm <= termination
        if (l2_norm <= l2_termination)
            break;
    }

    delete[] D;
    delete[] invD;
    delete[] R;
    delete[] local_y;
    delete[] res;

    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);

    return;
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