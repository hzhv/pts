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

std::vector<std::vector<int>> levelScheduling(const CSRMatrix &A);

void parallelTriangularSolve(
    const CSRMatrix &L, const std::vector<double> &b, std::vector<double> &x, 
    const std::vector<std::vector<int>> &levels, int rank, int numProcs
);

void serialTriangularSolve(
    const CSRMatrix &L, const std::vector<double> &b, std::vector<double> &x
);

#endif