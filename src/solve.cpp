#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <numeric>
#include "common.h"

using namespace std;

//--------------------------------------------------
// broadcast Levels
//--------------------------------------------------
void broadcastLevels(vector<vector<int>> &levels, int rank) {
    int numLevels = 0;
    if (rank == 0) {
        numLevels = levels.size();
    }
    MPI_Bcast(&numLevels, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        levels.resize(numLevels);
    }
    for (int i = 0; i < numLevels; i++) {
        int levelSize = 0;
        if (rank == 0) {
            levelSize = levels[i].size();
        }
        MPI_Bcast(&levelSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            levels[i].resize(levelSize);
        }
        if (levelSize > 0) {
            MPI_Bcast(levels[i].data(), levelSize, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }
}

//---------------------------------------------------------------
// BlockRowPartition: divide the global CSR matrix into blocks
//---------------------------------------------------------------
void blockRowPartition(const CSRMatrix& A_csr, vector<LocalCSR>& localParts, int numProcs)
{
    int n = A_csr.n;
    int baseSize = n / numProcs;
    int remainder = n % numProcs;
    auto rowBeginOf = [&](int i) {
        return i < remainder ? i*(baseSize+1) : i*baseSize + remainder;
    };
    localParts.resize(numProcs);
    for (int r = 0; r < numProcs; r++) {
        LocalCSR part;
        part.start_global_row = rowBeginOf(r);
        int rowEnd = (r < remainder) ? part.start_global_row + (baseSize + 1)
                                     : part.start_global_row + baseSize;
        part.local_n = rowEnd - part.start_global_row;
        part.row_ptr.resize(part.local_n + 1, 0);
        int nnz_local = 0;
        for (int i = 0; i < part.local_n; i++) {
            int globalRow = part.start_global_row + i;
            int rowCount = A_csr.row_ptr[globalRow+1] - A_csr.row_ptr[globalRow];
            part.row_ptr[i+1] = rowCount;
            nnz_local += rowCount;
        }
        for (int i = 0; i < part.local_n; i++) {
            part.row_ptr[i+1] += part.row_ptr[i];
        }
        part.col_id.resize(nnz_local);
        part.val.resize(nnz_local);
        for (int i = 0; i < part.local_n; i++) {
            int globalRow = part.start_global_row + i;
            int rowStart = A_csr.row_ptr[globalRow];
            int rowEnd2  = A_csr.row_ptr[globalRow+1];
            int localOffset = part.row_ptr[i];
            for (int idx = rowStart; idx < rowEnd2; idx++) {
                int locPos = localOffset + (idx - rowStart);
                part.col_id[locPos] = A_csr.col_id[idx];
                part.val[locPos]    = A_csr.val[idx];
            }
        }
        localParts[r] = std::move(part);
    }
}

//--------------------------------------------------
// distributeBlocks: communicate by MPI_Scatterv, 分发全局 CSR 数据到各进程
// root (rank 0) has globalMatrix
//--------------------------------------------------
void distributeBlocks(const CSRMatrix& globalMatrix, LocalCSR& localMatrix, int rank, int numProcs) {
    int n = globalMatrix.n;
    vector<int> rowsPerProc(numProcs, 0);
    vector<int> displsRows(numProcs, 0);
    int base = n / numProcs;
    int rem = n % numProcs;
    for (int i = 0; i < numProcs; i++){
         rowsPerProc[i] = base + (i < rem ? 1 : 0);
         if (i > 0)
             displsRows[i] = displsRows[i-1] + rowsPerProc[i-1];
    }
    vector<int> sendcounts_row(numProcs), displs_row(numProcs);
    for (int i = 0; i < numProcs; i++){
        sendcounts_row[i] = rowsPerProc[i] + 1;
        displs_row[i] = displsRows[i];
    }
    vector<int> sendcounts_data(numProcs), displs_data(numProcs);
    for (int i = 0; i < numProcs; i++){
         int start = displsRows[i];
         int end = start + rowsPerProc[i];
         sendcounts_data[i] = globalMatrix.row_ptr[end] - globalMatrix.row_ptr[start];
         displs_data[i] = globalMatrix.row_ptr[start];
    }
    int local_row_count = rowsPerProc[rank];
    int local_nnz = sendcounts_data[rank];
    localMatrix.local_n = local_row_count;
    localMatrix.start_global_row = displsRows[rank];
    localMatrix.row_ptr.resize(local_row_count + 1);
    localMatrix.col_id.resize(local_nnz);
    localMatrix.val.resize(local_nnz);
    // 仅 root 使用全局数据，其他进程传入 nullptr
    int* sendbuf_row = (rank == 0) ? const_cast<int*>(globalMatrix.row_ptr.data()) : nullptr;
    MPI_Scatterv(sendbuf_row, sendcounts_row.data(), displs_row.data(), MPI_INT,
                 localMatrix.row_ptr.data(), local_row_count + 1, MPI_INT,
                 0, MPI_COMM_WORLD);
    int* sendbuf_col = (rank == 0) ? const_cast<int*>(globalMatrix.col_id.data()) : nullptr;
    MPI_Scatterv(sendbuf_col, sendcounts_data.data(), displs_data.data(), MPI_INT,
                 localMatrix.col_id.data(), local_nnz, MPI_INT,
                 0, MPI_COMM_WORLD);
    double* sendbuf_val = (rank == 0) ? const_cast<double*>(globalMatrix.val.data()) : nullptr;
    MPI_Scatterv(sendbuf_val, sendcounts_data.data(), displs_data.data(), MPI_DOUBLE,
                 localMatrix.val.data(), local_nnz, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    // 调整本地 row_ptr 使其以 0 为起始值
    int base_val = localMatrix.row_ptr[0];
    for (int i = 0; i < (int)localMatrix.row_ptr.size(); i++){
         localMatrix.row_ptr[i] -= base_val;
    }
}

//--------------------------------------------------
// Block Based Parallel Triangular Solver（改用 MPI_Allgatherv 同步局部解）
//--------------------------------------------------
void parallelTriangularSolve_block(const LocalCSR& myPart, const vector<double>& b, 
    vector<double>& globalX,  
    const vector<vector<int>> &levels, int rank, int numProcs)
{
    int localStart = myPart.start_global_row;
    int localEnd   = localStart + myPart.local_n;
    vector<double> localX(myPart.local_n, 0.0);

    // 对于 MPI_Allgatherv，先收集每个进程的局部解长度
    vector<int> recvcounts(numProcs, 0);
    MPI_Allgather(&myPart.local_n, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    vector<int> displs(numProcs, 0);
    for (int i = 1; i < numProcs; i++) {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    // 遍历每一层（levels 存储全局行号）
    for (size_t lev = 0; lev < levels.size(); lev++) {
        for (int row : levels[lev]) {
            if (row >= localStart && row < localEnd) {
                int localRow = row - localStart;
                double sum = 0.0;
                int rowBegin = myPart.row_ptr[localRow];
                int rowEnd2  = myPart.row_ptr[localRow+1];
                for (int idx = rowBegin; idx < rowEnd2; idx++) {
                    int j = myPart.col_id[idx];
                    if (j < row) {
                        sum += myPart.val[idx] * globalX[j];
                    }
                }
                double diag = 0.0;
                for (int idx = rowBegin; idx < rowEnd2; idx++) {
                    if (myPart.col_id[idx] == row) {
                        diag = myPart.val[idx];
                        break;
                    }
                }
                if (diag == 0.0) {
                    throw runtime_error("Zero diagonal at row " + to_string(row));
                }
                localX[localRow] = (b[row] - sum) / diag;
            }
        }
        // 将本地计算结果写回 globalX 对应部分
        for (int i = localStart; i < localEnd; i++) {
            globalX[i] = localX[i - localStart];
        }
        // 同步所有进程的局部解到全局解向量 globalX 使用 MPI_Allgatherv
        MPI_Allgatherv(globalX.data() + localStart, myPart.local_n, MPI_DOUBLE,
                       globalX.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                       MPI_COMM_WORLD);
    }
}

// void distributeBlocks(const std::vector<LocalCSR>& localParts, 
//                       LocalCSR& myPart, 
//                       int rank, int numProcs)
// {
//     // rank 0 给自己直接赋值
//     if (rank == 0) {
//     myPart = localParts[0];
//     }

//     // 其余进程：rank 0 用 MPI_Send 发, rank r 用 MPI_Recv 收
//     // row_ptr
//     if (rank == 0) {
//         for (int r = 1; r < numProcs; r++) {
//             const auto& p = localParts[r];
//             int local_n = p.local_n;
//             int row_ptr_size = p.row_ptr.size();
//             int col_id_size  = p.col_id.size();
//             int val_size     = p.val.size();

//             // 1) 先发送 local_n, start_global_row
//             MPI_Send(&local_n, 1, MPI_INT, r, 101, MPI_COMM_WORLD);
//             MPI_Send(&p.start_global_row, 1, MPI_INT, r, 102, MPI_COMM_WORLD);
//             // 2) row_ptr
//             MPI_Send(&row_ptr_size, 1, MPI_INT, r, 103, MPI_COMM_WORLD);
//             MPI_Send(p.row_ptr.data(), row_ptr_size, MPI_INT, r, 104, MPI_COMM_WORLD);
//             // 3) col_id
//             MPI_Send(&col_id_size, 1, MPI_INT, r, 105, MPI_COMM_WORLD);
//             MPI_Send(p.col_id.data(), col_id_size, MPI_INT, r, 106, MPI_COMM_WORLD);
//             // 4) val
//             MPI_Send(&val_size, 1, MPI_INT, r, 107, MPI_COMM_WORLD);
//             MPI_Send(p.val.data(), val_size, MPI_DOUBLE, r, 108, MPI_COMM_WORLD);
//         }
//     }
//     else {
//         // rank != 0: 接收
//         int local_n;
//         MPI_Recv(&local_n, 1, MPI_INT, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         MPI_Recv(&myPart.start_global_row, 1, MPI_INT, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         myPart.local_n = local_n;

//         int row_ptr_size, col_id_size, val_size;
//         // row_ptr
//         MPI_Recv(&row_ptr_size, 1, MPI_INT, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         myPart.row_ptr.resize(row_ptr_size);
//         MPI_Recv(myPart.row_ptr.data(), row_ptr_size, MPI_INT, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         // col_id
//         MPI_Recv(&col_id_size, 1, MPI_INT, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         myPart.col_id.resize(col_id_size);
//         MPI_Recv(myPart.col_id.data(), col_id_size, MPI_INT, 0, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         // val
//         MPI_Recv(&val_size, 1, MPI_INT, 0, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         myPart.val.resize(val_size);
//         MPI_Recv(myPart.val.data(), val_size, MPI_DOUBLE, 0, 108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }
// }


void parallelTriangularSolve(const CSRMatrix& L, const vector<double>& b, vector<double>& x, 
                               const vector<vector<int>>& levels, int rank, int numProcs) 
{
    int n = L.n;
    x.assign(n, 0.0);

    for (size_t lev = 0; lev < levels.size(); lev++) {
        for (int row : levels[lev]) {
            if ((row % numProcs) == rank) {
                double sum = 0.0;
                for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; idx++) {
                    int j = L.col_id[idx];
                    if (j < row) {
                        sum += L.val[idx] * x[j];
                    }
                }
                double diag = 0.0;
                for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; idx++) {
                    if (L.col_id[idx] == row) {
                        diag = L.val[idx];
                        break;
                    }
                }
                if (diag == 0.0) {
                    throw runtime_error("Zero diagonal encountered at row " + to_string(row));
                }
                x[row] = (b[row] - sum) / diag;
            }
            MPI_Bcast(&x[row], 1, MPI_DOUBLE, row % numProcs, MPI_COMM_WORLD);
        }
    }
}


// ------------------ The following is for serial solving and Testing only --------------------
void serialTriangularSolve(const CSRMatrix& L, const vector<double>& b, vector<double>& x) 
{
    int n = L.n;
    x.assign(n, 0.0);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        double diag = 0.0;
        for (int idx = L.row_ptr[i]; idx < L.row_ptr[i+1]; idx++) {
            int col = L.col_id[idx];
            if (col < i) {
                sum += L.val[idx] * x[col];
            } else if (col == i) {
                diag = L.val[idx];
            }
        }
        if (diag == 0.0) {
            throw runtime_error("Zero diagonal encountered at row " + to_string(i));
        }
        x[i] = (b[i] - sum) / diag;
    }
}
// ------------------------------------------------------------------------------------------
