#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "common.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <triplet_filename>\n";
        return -1;
    }

    MPI_Init(&argc, &argv);
    
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    CSRMatrix A_csr;
    int row_ptr_size = 0;
    int col_id_size = 0;
    int val_size = 0;
    std::vector<std::vector<int>> levels;

    if (rank == 0) {
        std::string CSR_file = argv[1];
        std::ifstream infile(CSR_file);
        if (!infile) {
            std::cerr << "Cannot open file: " << CSR_file << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        std::cout << "Loading A in CSR format..." << std::endl;
        A_csr = loadCSR(CSR_file);
        std::cout << "CSR Matrix loaded:" << std::endl;
        std::cout << "Number of rows: " << A_csr.n << std::endl;
        std::cout << "Nonzeros: " << A_csr.val.size() << std::endl;

        row_ptr_size = A_csr.row_ptr.size();
        col_id_size = A_csr.col_id.size();
        val_size = A_csr.val.size();
    }

    // Bcast A_csr
    MPI_Bcast(&A_csr.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&row_ptr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col_id_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&val_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        A_csr.row_ptr.resize(row_ptr_size);
        A_csr.col_id.resize(col_id_size);
        A_csr.val.resize(val_size);
    }

    MPI_Bcast(A_csr.row_ptr.data(), row_ptr_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_csr.col_id.data(), col_id_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_csr.val.data(), val_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Broadcast of CSRMatrix completed. Matrix n = " << A_csr.n << std::endl;
        std::cout << std::endl;
    }

    std::vector<double> b(A_csr.n, 0.0);
    std::vector<double> x(A_csr.n, 0.0);
    // Bcast x, b
    if (rank == 0) {
        std::fill(b.begin(), b.end(), 1.0);  // b 全部初始化为 1
        std::fill(x.begin(), x.end(), 0.0);    // x 初始化为 0
    }
    MPI_Bcast(b.data(), A_csr.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x.data(), A_csr.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    try {   
        levels = levelScheduling(A_csr); // Not finished, need to be bcast in the future
        if (rank == 0)
            std::cout << "There are " << levels.size() << " levels in total." << std::endl;

        std::cout << "Rank " << rank << " entering barrier before parallelTriangularSolve." << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "Rank " << rank << " leaving barrier before parallelTriangularSolve." << std::endl;
        std::cout << std::endl;

        double start_parallel = MPI_Wtime();
        parallelTriangularSolve(A_csr, b, x, levels, rank, numProcs);
        double end_parallel = MPI_Wtime();
        
        if (rank == 0) {
            std::cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << std::endl;
            std::cout << std::endl;

            std::ofstream outfile("Parallel_solution.txt");
            if (!outfile) {
                std::cerr << "Error: cannot open file for writing solution.\n";
                return -1;
            }
            for (int i = 0; i < A_csr.n; i++) {
                outfile << x[i] << "\n";
            }
            outfile.close();
            std::cout << "Solution x Saved.\n";
        }
    } 
    catch (const std::exception& e) {
        if (rank == 0)
            std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    MPI_Finalize();
    return 0;
}




// Hardcoded test case
// int main(int argc, char **argv) {
//     using namespace std;
//     MPI_Init(&argc, &argv);
    
//     int rank, numProcs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

//     try {
//         CSRMatrix L;
//         L.n = 8;
//         L.val = std::vector<double>(14, 1.0); // 
//         L.col_id = {
//             0,
//             1,
//             2,
//             0, 3,
//             0, 4,
//             1, 5,
//             2, 6,
//             3, 4, 7
//         };
//         L.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 15};
//         vector<double> b(L.n, 1.0); // b initialized to 1
//         vector<double> x(L.n, 0.0); // x initialized to 0
        
//         vector<vector<int>> levels = levelScheduling(L);
//         if (rank == 0) {
//             cout << "Level scheduling result:" << std::endl;
//             for (size_t lvl = 0; lvl < levels.size(); lvl++) {
//                 cout << "Level " << lvl << ": ";
//                 for (int row : levels[lvl])
//                     cout << row << " ";
//                 cout << std::endl;
//             }
//         }

//         MPI_Barrier(MPI_COMM_WORLD);
//         double start_parallel = MPI_Wtime();
//         parallelTriangularSolve(L, b, x, levels, rank, numProcs);
//         double end_parallel = MPI_Wtime();
        
//         if (rank == 0) {
//             cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << std::endl;
//             cout << "Parallel solution x:" << std::endl;
//             for (int i = 0; i < L.n; i++) {
//                 cout << x[i] << " ";
//             }
//             cout << std::endl;
//         }
        
//         // Serial solve for comparison (only on rank 0)
//         cout << '\n';
//         if (rank == 0) {
//             vector<double> x_serial(L.n, 0.0);
//             double start_serial = MPI_Wtime();
//             serialTriangularSolve(L, b, x_serial);
//             double end_serial = MPI_Wtime();
//             cout << "Serial solve time: " << (end_serial - start_serial) << " sec" << std::endl;
//             cout << "Serial solution x:" << std::endl;
//             for (int i = 0; i < L.n; i++) {
//                 cout << x_serial[i] << " ";
//             }
//             cout << std::endl;
//         }

//     } 
//     catch (const exception &e) {
//         if (rank == 0)
//             cerr << "Error: " << e.what() << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, -1);
//     }
    
//     MPI_Finalize();
//     return 0;
// }