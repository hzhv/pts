#include <mpi.h>
#include <vector>
#include <algorithm>
#include <numeric> 
#include <stdexcept>
#include "common.h"

struct RowRange 
{ 
    int first;
    int last; 
};

static std::vector<RowRange> makeBlockPartition(int P, int n) 
/**
 * @brief Partition the rows [0, n) into P contiguous blocks as evenly as possible.
 * e.g. if n = 10 and P = 3, the result will be:
 *      Block 0: [0, 4)
 *      Block 1: [4, 7)
 *      Block 2: [7, 10)
 * @param P The number of blocks (e.g., number of processes).
 * @param n The total number of elements to partition.
 * @return A vector of RowRange structs, each representing the range assigned to one block.
 */
{
    std::vector<RowRange> ranges(P);
    int base = n / P;
    int rem = n % P; 
    int start = 0;
    for (int b = 0; b < P; ++b) {
        int len = base + (b < rem ? 1 : 0);
        ranges[b] = { start, start + len };
        start += len;
    }
    return ranges;
}

void parallelTriangularSolve_block(
    const CSRMatrix& L,
    std::vector<double>& x, const std::vector<double>& b,
    std::vector<int>& level_ptr,
    std::vector<int>& level_rows,
    std::vector<int>& row_owner,
    std::vector<int>& dep_ptr,
    std::vector<int>& dep_rows,
    int rank, int P
) 
/**
 * @brief Parallel triangular solve using block partitioning and Allgatherv.
 * 
 */
{
    int n = L.n;
    int numLevels =(int)level_ptr.size() - 1;

    // auto ranges = makeBlockPartition(P, n);
    // auto owner = [&](int row)
    // {
    //     // Changed from row % P
    //     for (int rank = 0; rank < P; ++rank)
    //         if (row >= ranges[rank].first && row < ranges[rank].last)
    //             return rank;
    //     return 0; 

    //     // return row % P;
    // };

    // local_level_ptr[k]..[k+1]) 存放本 rank 在第 k 层的行在 local_level_rows 中的区间
    std::vector<int> local_level_ptr(numLevels+1, 0);
    std::vector<int> local_level_rows;
    // level_counts[k][r] = kth level's row count on rank r
    std::vector<std::vector<int>> level_counts(numLevels, std::vector<int>(P, 0)); // lvls x ranks
    for (int k = 0; k < numLevels; ++k) 
    {
        for (int pos = level_ptr[k]; pos < level_ptr[k+1]; ++pos) 
        {
            int row = level_rows[pos];
            int r = row_owner[row];
            level_counts[k][r] += 1;
            if (r == rank)
                local_level_rows.push_back(row);
        }
        local_level_ptr[k+1] = (int)local_level_rows.size();
    }

    // Each levels' recvbuf position for MPI_Allgatherv
    std::vector<std::vector<int>> level_displs(numLevels, std::vector<int>(P, 0));
    for (int k = 0; k < numLevels; ++k) 
    {
        int acc = 0;
        for (int r = 0; r < P; ++r) 
        {
            level_displs[k][r] = acc;
            acc += level_counts[k][r];
        }
    }

// ========================== New Test ===============================
// 
    // --- 2) 一次性构造 send_rows_ptr & send_rows ---
    //    直接重用 local_level_ptr & local_level_rows
    // const auto& send_rows_ptr = local_level_ptr;
    // const auto& send_rows     = local_level_rows;

    // // --- 3) 一次性构造 recv_perm_ptr & recv_perm ---
    // std::vector<int> recv_perm_ptr(numLevels+1, 0);
    // std::vector<int> recv_perm;
    // recv_perm.reserve(level_rows.size());
    // for (int k = 0; k < numLevels; ++k) {
    //     recv_perm_ptr[k] = (int)recv_perm.size();
    //     // 按 rank-块顺序，把本层行号依次写入 recv_perm
    //     for (int r = 0; r < P; ++r) {
    //         for (int pos = level_ptr[k]; pos < level_ptr[k+1]; ++pos) {
    //             int row = level_rows[pos];
    //             if (row_owner[row] == r) {
    //                 recv_perm.push_back(row);
    //             }
    //         }
    //     }
    // }
    // recv_perm_ptr[numLevels] = (int)recv_perm.size();

    // // 1d 
    // std::vector<int> flatCounts(numLevels * P), flatDispls(numLevels * P);
    // for (int k = 0; k < numLevels; ++k) {
    //     // counts
    //     for (int r = 0; r < P; ++r)
    //         flatCounts[k*P + r] = level_counts[k][r];
    //     // displs: prefix sum in each level
    //     int acc = 0;
    //     for (int r = 0; r < P; ++r) {
    //         flatDispls[k*P + r] = acc;
    //         acc += level_counts[k][r];
    //     }
    // }

    // // --- 5) 预分配 sendbuf/recvbuf 最大容量，避免循环 malloc ---
    // int maxSend = local_level_ptr.back();           // 本 rank 总行数
    // int maxRecv = level_ptr.back() - level_ptr[0];  // 最多整层行数 = n
    // std::vector<double> sendbuf, recvbuf;
    // sendbuf .reserve(maxSend);
    // recvbuf .reserve(maxRecv);

    // // --- 6) 阻塞 Allgatherv + 线性写回 的求解主循环 ---
    // x.assign(n, 0.0);
    // for (int k = 0; k < numLevels; ++k) {
    //     int sbeg    = send_rows_ptr[k];
    //     int sendCnt = send_rows_ptr[k+1] - sbeg;
    //     sendbuf.resize(sendCnt);
    //     for (int j = 0; j < sendCnt; ++j) {
    //         int row = send_rows[sbeg + j];
    //         double sum = 0.0, diag = 0.0;
    //         for (int p = L.row_ptr[row]; p < L.row_ptr[row+1]; ++p) {
    //             int c = L.col_id[p];
    //             double v = L.val[p];
    //             if (c == row) diag = v;
    //             else          sum += v * x[c];
    //         }
    //         sendbuf[j] = (b[row] - sum) / diag;
    //     }

    //     int Lb        = level_ptr[k];
    //     int Le        = level_ptr[k+1];
    //     int levelSize = Le - Lb;
    //     recvbuf.resize(levelSize);
    //     MPI_Allgatherv(
    //         sendbuf.data(), sendCnt, MPI_DOUBLE,
    //         recvbuf.data(),
    //         &flatCounts[k*P],
    //         &flatDispls[k*P],
    //         MPI_DOUBLE,
    //         MPI_COMM_WORLD
    //     );

    //     int base = recv_perm_ptr[k];
    //     for (int i = 0; i < levelSize; ++i) {
    //         x[ recv_perm[base + i] ] = recvbuf[i];
    //     }
    // }
// ===================================================================================
    // // Main loop: Each level do a Allgatherv
    for (int k = 0; k < numLevels; ++k) {
        int Lb = level_ptr[k], Le = level_ptr[k+1];
        int levelSize = Le - Lb;  // current level size

        int lb = local_level_ptr[k], le = local_level_ptr[k+1];
        int localCount = le - lb; // equivalent to level_counts[k][rank]

        // Calc local x
        std::vector<double> sendbuf(localCount);
        for (int idx = lb; idx < le; ++idx) 
        {
            int row = local_level_rows[idx];
            double sum = 0.0;
            int diag_idx = L.row_ptr[row+1] - 1;
            double diag = L.val[diag_idx];
            for (int p = L.row_ptr[row]; p < L.row_ptr[row+1]; ++p)
            {
                int col = L.col_id[p];
                sum  += L.val[p] * x[col];
            }
            sendbuf[idx - lb] = (b[row] - sum) / diag;
        }

        std::vector<double> recvbuf(levelSize);
        // // ----- Blocked
        MPI_Allgatherv(
            sendbuf.data(), localCount, MPI_DOUBLE,
            recvbuf.data(),
            level_counts[k].data(),
            level_displs[k].data(),
            MPI_DOUBLE,
            MPI_COMM_WORLD
        );

        int pos = 0;
        for (int r = 0; r < P; ++r) {
            // 针对本层，找出所有 owner(row)==r 的行号，顺序与 level_rows 相同
            for (int idx = Lb; idx < Le; ++idx) {
                int row = level_rows[idx];
                if (row_owner[row] == r) {          
                    x[row] = recvbuf[pos++];    
                }
            }
        }
    }

}


