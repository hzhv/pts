#pragma once

#include <vector>

struct Message {
    int row;       // 已计算行的全局行号
    double value;  // 对应行计算后的解 x[row]
};

struct CSRMatrix 
{
    int n;                    
    std::vector<int> row_ptr;      
    std::vector<int> col_id;      
    std::vector<double> val;       
};

struct Triplet {
    int row;
    int col;
    double val;
};

CSRMatrix loadCSR(std::string& filename);

CSRMatrix loadCSRFromTripletFile(const std::string& filename);

std::vector<std::vector<int>> levelScheduling
(
    const CSRMatrix& A, 
    std::vector<std::vector<int>>& dependents,
    std::vector<std::vector<int>>& dependencies
);

void levelScheduling_plain
(
    const CSRMatrix& A,
    std::vector<int>& level_ptr,
    std::vector<int>& level_rows,
    std::vector<int>& dep_ptr,
    std::vector<int>& dep_rows
);

void buildSchedules(
    const CSRMatrix& A,
    const std::vector<int>& row_owner,   
    int rank,
    int P,

    // ---- 输出 1：层信息
    std::vector<int>& level_ptr,
    std::vector<int>& level_rows,

    // ---- 输出 2：Allgatherv 计数 + 偏移 (flattened)
    std::vector<int>& level_counts_flat,      // (#levels)*P
    std::vector<int>& level_displs_flat,

    // ---- 输出 3：本 rank 需要发送的行
    std::vector<int>& send_rows_ptr,          // (#levels)+1
    std::vector<int>& send_rows,

    // ---- 输出 4：recvbuf -> x 的置换
    std::vector<int>& recv_perm_ptr,          // (#levels)+1
    std::vector<int>& recv_perm
);

void parallelTriangularSolve
(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x, 
    const std::vector<std::vector<int>>& levels, int rank, int numProcs
);

void parallelTriangularSolve_p2p
(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x, 
    const std::vector<std::vector<int>>& dependents,
    const std::vector<std::vector<int>>& dependencies,
    const std::vector<std::vector<int>>& levels, int rank, int numProcs
);

void parallelTriangularSolve_block
(
    const CSRMatrix& L,
    std::vector<double>& x, const std::vector<double>& b,
    std::vector<int>& level_ptr,
    std::vector<int>& level_rows,
    std::vector<int>& row_owner,
    std::vector<int>& dep_ptr,
    std::vector<int>& dep_rows,
    int rank, int P
    // MPI_Comm comm = MPI_COMM_WORLD
);

void parallelTriangularSolve_fast
(
    const CSRMatrix& L,
    const std::vector<int>& level_ptr,
    const std::vector<int>& level_rows,
    const std::vector<int>& level_counts_flat,  
    const std::vector<int>& level_displs_flat,  
    const std::vector<int>& send_rows_ptr,
    const std::vector<int>& send_rows,          
    const std::vector<int>& recv_perm_ptr,
    const std::vector<int>& recv_perm,
    const std::vector<double>& b,
    std::vector<double>& x,
    int rank, int P
);

void serialTriangularSolve
(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x
);
