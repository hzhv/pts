#ifndef COMMON_H
#define COMMON_H

#include <vector>

struct CSRMatrix 
{
    int n;                    
    std::vector<int> row_ptr;      
    std::vector<int> col_id;      
    std::vector<double> val;       
};

struct LocalCSR 
{
    int local_n;             // 本地拥有的行数
    int start_global_row;    // 该块在全局的起始行号
    std::vector<int> row_ptr;
    std::vector<int> col_id;
    std::vector<double> val;
};

struct Triplet {
    int row;
    int col;
    double val;
};

std::vector<std::vector<int>> levelScheduling(const CSRMatrix& A);

CSRMatrix loadCSRFromTripletFile(const std::string& filename);
CSRMatrix loadCSR(std::string& filename);

void parallelTriangularSolve(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x, 
    const std::vector<std::vector<int>>& levels, int rank, int numProcs
);

void serialTriangularSolve(
    const CSRMatrix& L, const std::vector<double>& b, std::vector<double>& x
);

// void distributeBlocks(const std::vector<LocalCSR>& localParts, 
//     LocalCSR& myPart, 
//     int rank, int numProcs);
void distributeBlocks(const CSRMatrix& globalMatrix, LocalCSR& localMatrix, 
    int rank, int numProcs);

void parallelTriangularSolve_block(const LocalCSR& myPart, const std::vector<double>& b, 
        std::vector<double>& globalX,
        const std::vector<std::vector<int>> &levels, int rank, int numProcs);

void broadcastLevels(std::vector<std::vector<int>>& levels, int rank);

void blockRowPartition(const CSRMatrix& A_csr, std::vector<LocalCSR>& localParts, int numProcs);

#endif