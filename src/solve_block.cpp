#include <mpi.h>
#include <vector>
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

    auto ranges = makeBlockPartition(P, n);
    auto owner = [&](int row)
    {
        // Changed from row % P
        for (int rank = 0; rank < P; ++rank)
            if (row >= ranges[rank].first && row < ranges[rank].last)
                return rank;
        return 0; 

        // return row % P;
    };

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
            int r = owner(row);
            level_counts[k][r] += 1;
            if (r == rank)
                local_level_rows.push_back(row);
        }
        local_level_ptr[k+1] = (int)local_level_rows.size();
    }

    //预计算每层的 offset for Allgatherv
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

    // Main loop: Each level do a Allgatherv
    for (int k = 0; k < numLevels; ++k) {
        int Lb = level_ptr[k], Le = level_ptr[k+1];
        int levelSize = Le - Lb;  // current level size

        // 本 rank 在本层的本地行区间
        int lb = local_level_ptr[k], le = local_level_ptr[k+1];
        int localCount = le - lb;

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
        // MPI_Allgatherv(
        //     sendbuf.data(), localCount, MPI_DOUBLE,
        //     recvbuf.data(),
        //     level_counts[k].data(),
        //     level_displs[k].data(),
        //     MPI_DOUBLE,
        //     MPI_COMM_WORLD
        // );

        // ----- Non-blocked
        MPI_Request req;
        MPI_Iallgatherv(
            sendbuf.data(), localCount, MPI_DOUBLE,
            recvbuf.data(),
            level_counts[k].data(),
            level_displs[k].data(),
            MPI_DOUBLE,
            MPI_COMM_WORLD,
            &req
        );
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        // for (int i = 0; i < levelSize; ++i) {
        //     int row = level_rows[Lb + i];
        //     x[row] = recvbuf[i];
        // }

        int pos = 0;
        for (int r = 0; r < P; ++r) {
            // 针对本层，找出所有 owner(row)==r 的行号，顺序与 level_rows 相同
            for (int idx = Lb; idx < Le; ++idx) {
                int row = level_rows[idx];
                if (owner(row) == r) {          
                    x[row] = recvbuf[pos++];    
                }
            }
        }
    }
}


