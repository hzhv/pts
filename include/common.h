#ifndef COMMON_H
#define COMMON_H

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
    std::vector<int>& dep_ptr,
    std::vector<int>& dep_rows,
    int rank, int P
    // MPI_Comm comm = MPI_COMM_WORLD
);

void serialTriangularSolve
(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x
);

// void distributeBlocks(const std::vector<LocalCSR>& localParts, 
//     LocalCSR& myPart, 
//     int rank, int numProcs);
// struct LocalCSR 
// {
//     int local_n;             // 本地拥有的行数
//     int start_global_row;    // 该块在全局的起始行号
//     std::vector<int> row_ptr;
//     std::vector<int> col_id;
//     std::vector<double> val;
// };

// void distributeBlocks(const CSRMatrix& globalMatrix, LocalCSR& localMatrix, 
//     int rank, int numProcs);

// void parallelTriangularSolve_block(const LocalCSR& myPart, const std::vector<double>& b, 
//         std::vector<double>& globalX,
//         const std::vector<std::vector<int>> &levels, int rank, int numProcs);

// void broadcastLevels(std::vector<std::vector<int>>& levels, int rank);

// void blockRowPartition(const CSRMatrix& A_csr, std::vector<LocalCSR>& localParts, int numProcs);

#endif