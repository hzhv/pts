#include <mpi.h>
#include <vector>
#include <numeric> 
#include <stdexcept>
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
    int rank, int P)
{
    const int numLevels = (int)level_ptr.size() - 1;
    std::vector<double> sendbuf, recvbuf;
    for (int k = 0; k < numLevels; ++k) {

        // ----- local computation
        int sbeg = send_rows_ptr[k];
        int sendCnt = send_rows_ptr[k+1] - sbeg;
        sendbuf.resize(sendCnt);
        for (int j = 0; j < sendCnt; ++j) {
            int row = send_rows[sbeg + j];
            double sum = 0.0, diag = 0.0;
            for (int p = L.row_ptr[row]; p < L.row_ptr[row+1]; ++p) {
                int col = L.col_id[p];
                double v = L.val[p];
                if (col == row) diag = v;
                else            sum += v * x[col];
            }
            sendbuf[j] = (b[row] - sum) / diag;
        }

        int rbeg = level_ptr[k], rend = level_ptr[k+1];
        int levelSize = rend - rbeg;
        // // ----- Blocked
        // recvbuf.resize(levelSize);
        // MPI_Allgatherv(sendbuf.data(), sendCnt, MPI_DOUBLE,
        //                recvbuf.data(),
        //                &level_counts_flat[k*P],
        //                &level_displs_flat[k*P],
        //                MPI_DOUBLE, MPI_COMM_WORLD);
        
        // ----- Non-blocked
        recvbuf.resize(levelSize);
        MPI_Request req;
        MPI_Iallgatherv(
            sendbuf.data(), sendCnt, MPI_DOUBLE,
            recvbuf.data(),
            &level_counts_flat[k*P],
            &level_displs_flat[k*P],
            MPI_DOUBLE, MPI_COMM_WORLD,
            &req
        );
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        // ----- write-back  (no branch)
        int permBase = recv_perm_ptr[k];
        for (int i = 0; i < levelSize; ++i)
            x[ recv_perm[permBase + i] ] = recvbuf[i];
    }
}
// --mca coll_sync_barrier false
