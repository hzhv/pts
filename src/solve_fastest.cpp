#include <mpi.h>
#include <vector>
#include <algorithm>
#include <numeric> 
#include <stdexcept>
#include <iostream>   
#include "common.h"

void parallelTriangularSolve_fast(
    const CSRMatrix& L,
    const std::vector<int>& level_ptr,
    const std::vector<int>& level_rows,
    const std::vector<int>& level_counts_flat,  // (#levels)*P
    const std::vector<int>& level_displs_flat,  // (#levels)*P
    const std::vector<int>& send_rows_ptr,
    const std::vector<int>& send_rows,          // 本 rank 行
    const std::vector<int>& recv_perm_ptr,
    const std::vector<int>& recv_perm,
    const std::vector<double>& b,
    std::vector<double>& x,
    double& comm_time,
    int rank, int P)
{
    const int numLevels = (int)level_ptr.size() - 1;
    std::vector<double> sendbuf, recvbuf;

    int maxSend = *std::max_element(send_rows_ptr.begin()+1, send_rows_ptr.end());
    int maxRecv = level_ptr.back();  // worst case: n

    sendbuf.reserve(maxSend);
    recvbuf.reserve(maxRecv);
    for (int k = 0; k < numLevels; ++k) {
        int sbeg     = send_rows_ptr[k];
        int sendCnt  = send_rows_ptr[k+1] - sbeg;
        sendbuf.resize(sendCnt);
        for (int j = 0; j < sendCnt; ++j) {
            int row = send_rows[sbeg + j];
            double sum = 0.0;
            int diag_idx = L.row_ptr[row+1] - 1;
            double diag = L.val[diag_idx];
            for (int p = L.row_ptr[row]; p < L.row_ptr[row+1]; ++p) {
                sum += L.val[p] * x[L.col_id[p]];
            }
            sendbuf[j] = (b[row] - sum) / diag;
        }

        int rbeg      = level_ptr[k];
        int rend      = level_ptr[k+1];
        int levelSize = rend - rbeg;
        recvbuf.resize(levelSize);

        // double t0 = MPI_Wtime();
        MPI_Allgatherv(sendbuf.data(), sendCnt, MPI_DOUBLE,
                       recvbuf.data(),
                       &level_counts_flat[k*P],
                       &level_displs_flat[k*P],
                       MPI_DOUBLE, MPI_COMM_WORLD);
        // comm_time += MPI_Wtime() - t0;

        int permBase = recv_perm_ptr[k];
        for (int i = 0; i < levelSize; ++i) {
            x[ recv_perm[permBase + i] ] = recvbuf[i];
        }
    }
}

// --mca coll_sync_barrier false
