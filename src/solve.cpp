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

void parallelTriangularSolve(const CSRMatrix& L, const vector<double>& b, vector<double>& x, 
                               const vector<vector<int>>& levels, int rank, int numProcs) 
{
    for (size_t lev = 0; lev < levels.size(); ++lev) {
        for (int row : levels[lev]) {
            if ((row % numProcs) == rank) {
                double sum = 0.0;
                int diag_idx = L.row_ptr[row+1] - 1;
                double diag = L.val[diag_idx];
                for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; ++idx) {
                    int col = L.col_id[idx];
                    sum += L.val[idx] * x[col];
                }
                // double diag = 0.0;
                // for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; idx++) {
                //     if (L.col_id[idx] == row) {
                //         diag = L.val[idx];
                //         break;
                //     }
                // }
                // if (diag == 0.0) {
                //     throw runtime_error("Zero diagonal encountered at row " + to_string(row));
                // }
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
// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <csr_filename>\n";
//         return -1;
//     }
//     MPI_Init(&argc, &argv);
//     int rank, numProcs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

//     // 1) Rank 0 load global L in CSR format
//     CSRMatrix A_csr;
//     if (rank == 0) {
//         std::string CSR_file = argv[1];
//         std::ifstream infile(CSR_file);
//         if (!infile) {
//             std::cerr << "Cannot open file: " << CSR_file << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, -1);
//         }
//         std::cout << "Loading A in CSR format..." << std::endl;
//         A_csr = loadCSR(CSR_file);
//         std::cout << "CSR Matrix loaded: n = " << A_csr.n 
//              << ", nnz = " << A_csr.val.size() << std::endl;
//     }

//     // 2) bcast # of rows, row_ptr
//     MPI_Bcast(&A_csr.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     int row_ptr_size = 0;
//     if (rank == 0) {
//         row_ptr_size = A_csr.row_ptr.size();
//     }
//     MPI_Bcast(&row_ptr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     if (rank != 0) {
//         A_csr.row_ptr.resize(row_ptr_size);
//     }
//     MPI_Bcast(A_csr.row_ptr.data(), row_ptr_size, MPI_INT, 0, MPI_COMM_WORLD);
//     // Note：other data col_id, val did not bcast here

//     // 3) Patition: root computes the partition info, then distributes to all processes
//     LocalCSR myPart;
//     // 这里采用 MPI_Scatterv 分发全局 CSR 数据（仅使用 A_csr.row_ptr, col_id, val 在 root 上）
//     distributeBlocks(A_csr, myPart, rank, numProcs);

//     // 4) init b , and globalX
//     std::vector<double> b;
//     int nGlobal;
//     if (rank == 0) {
//         nGlobal = A_csr.n;
//         b.resize(nGlobal, 1.0);
//     }
//     MPI_Bcast(&nGlobal, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     if (rank != 0) b.resize(nGlobal);
//     MPI_Bcast(b.data(), nGlobal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     std::vector<double> globalX(nGlobal, 0.0);

//     // 5) Level scheduling：only compute on rank 0, then bcast
//     std::vector<std::vector<int>> levels;
//     if (rank == 0) {
//         levels = levelScheduling(A_csr);
//         std::cout << "Level scheduling computed: " << levels.size() << " levels." << std::endl;
//     }
//     broadcastLevels(levels, rank);

//     std::cout << "Rank " << rank << " entering parallelTriangularSolve_block." << std::endl;
//     MPI_Barrier(MPI_COMM_WORLD);
//     std::cout << "Rank " << rank << " leaving parallelTriangularSolve_block." << std::endl;
//     double tstart = MPI_Wtime();
//     parallelTriangularSolve_block(myPart, b, globalX, levels, rank, numProcs);
//     double tend = MPI_Wtime();

//     if (rank == 0) {
//         std::cout << "Parallel solve time: " << (tend - tstart) << " s" << std::endl;
//         std::cout << std::endl;
//     }
//     MPI_Finalize();
//     return 0;
// }