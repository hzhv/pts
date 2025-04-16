#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "common.h"

// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <triplet_filename>\n";
//         return -1;
//     }

//     MPI_Init(&argc, &argv);
    
//     int rank, numProcs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

//     CSRMatrix A_csr;
//     int row_ptr_size = 0;
//     int col_id_size = 0;
//     int val_size = 0;
//     std::vector<std::vector<int>> levels;

//     if (rank == 0) {
//         std::string CSR_file = argv[1];
//         std::ifstream infile(CSR_file);
//         if (!infile) {
//             std::cerr << "Cannot open file: " << CSR_file << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, -1);
//         }
//         std::cout << "Loading A in CSR format..." << std::endl;
//         A_csr = loadCSR(CSR_file);
//         std::cout << "CSR Matrix loaded:" << std::endl;
//         std::cout << "Number of rows: " << A_csr.n << std::endl;
//         std::cout << "Nonzeros: " << A_csr.val.size() << std::endl;

//         row_ptr_size = A_csr.row_ptr.size();
//         col_id_size = A_csr.col_id.size();
//         val_size = A_csr.val.size();
//     }

//     // Bcast A_csr
//     MPI_Bcast(&A_csr.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&row_ptr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&col_id_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&val_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

//     if (rank != 0) {
//         A_csr.row_ptr.resize(row_ptr_size);
//         A_csr.col_id.resize(col_id_size);
//         A_csr.val.resize(val_size);
//     }

//     MPI_Bcast(A_csr.row_ptr.data(), row_ptr_size, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(A_csr.col_id.data(), col_id_size, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(A_csr.val.data(), val_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     std::vector<double> b(A_csr.n, 0.0);
//     std::vector<double> x(A_csr.n, 0.0);
//     // Bcast x, b
//     if (rank == 0) {
//         std::fill(b.begin(), b.end(), 1.0);    
//         std::fill(x.begin(), x.end(), 0.0);    
//     }
//     MPI_Bcast(b.data(), A_csr.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     MPI_Bcast(x.data(), A_csr.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     if (rank == 0) {
//         std::cout << "Broadcast of CSRMatrix completed. Matrix n = " << A_csr.n << std::endl;
//         std::cout << std::endl;
//     }

//     std::vector<std::vector<int>> dependents;
//     std::vector<std::vector<int>> dependencies;
//     std::vector<int> serialized_levels;
//     std::vector<int> serialized_dependents;
//     std::vector<int> serialized_dependencies;

//     if (rank == 0) {
//         levels = levelScheduling(A_csr, dependents, dependencies);
//         // Serialize levels, dependents, dependencies into a 1D array
//         for (const auto& level : levels) {
//             serialized_levels.push_back(level.size());
//             serialized_levels.insert(serialized_levels.end(), level.begin(), level.end());
//         }
//         for (int i = 0; i < A_csr.n; i++) {
//             serialized_dependents.push_back(dependents[i].size());
//             serialized_dependents.insert(serialized_dependents.end(), dependents[i].begin(), dependents[i].end());
//         }
//         for (int i = 0; i < A_csr.n; i++) {
//             serialized_dependencies.push_back(dependencies[i].size());
//             serialized_dependencies.insert(serialized_dependencies.end(), dependencies[i].begin(), dependencies[i].end());
//         }
//     }

//     // Bcast serialized dependents, dependencies, levels
//     int serialized_dependents_size = serialized_dependents.size();
//     MPI_Bcast(&serialized_dependents_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     if (rank != 0) {
//         serialized_dependents.resize(serialized_dependents_size);
//     }
//     MPI_Bcast(serialized_dependents.data(), serialized_dependents_size, MPI_INT, 0, MPI_COMM_WORLD);

//     int serialized_dependencies_size = serialized_dependencies.size();
//     MPI_Bcast(&serialized_dependencies_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     if (rank != 0) {
//         serialized_dependencies.resize(serialized_dependencies_size);
//     }
//     MPI_Bcast(serialized_dependencies.data(), serialized_dependencies_size, MPI_INT, 0, MPI_COMM_WORLD);

//     int serialized_size = serialized_levels.size();
//     MPI_Bcast(&serialized_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     if (rank != 0) {
//         serialized_levels.resize(serialized_size);
//     }
//     MPI_Bcast(serialized_levels.data(), serialized_size, MPI_INT, 0, MPI_COMM_WORLD);

//     // Deserialize dependents, dependencies, levels on all ranks
//     if (rank != 0) {
//         dependents.clear();
//         dependents.resize(A_csr.n);
//         for (int i = 0, idx = 0; i < A_csr.n; ++i) {
//             int count = serialized_dependents[idx++];
//             dependents[i].assign(serialized_dependents.begin() + idx, serialized_dependents.begin() + idx + count);
//             idx += count;
//         }

//         dependencies.clear();
//         dependencies.resize(A_csr.n);
//         for (int i = 0, idx = 0; i < A_csr.n; i++) {
//             int count = serialized_dependencies[idx++];
//             dependencies[i].assign(serialized_dependencies.begin() + idx, serialized_dependencies.begin() + idx + count);
//             idx += count;
//         }

//         levels.clear();
//         for (size_t i = 0; i < serialized_levels.size();) {
//             int level_size = serialized_levels[i++];
//             std::vector<int> level(serialized_levels.begin() + i, serialized_levels.begin() + i + level_size);
//             levels.push_back(level);
//             i += level_size;
//         }
//     }
//     if (rank == 0) {
//         std::cout << "Broadcast of dependents, dependencies, and levels completed." << std::endl;
//         std::cout << "There are " << levels.size() << " levels in total." << std::endl;
//         std::cout << std::endl;
//     }


//     try {
//         double start_parallel = MPI_Wtime();
//         // for (int i = 0; i < 100; ++i)
//         // parallelTriangularSolve(
//         //     A_csr, b, x,
//         //     levels, rank, numProcs);
//         parallelTriangularSolve_p2p(
//             A_csr, b, x,
//             dependents, 
//             dependencies,
//             levels, rank, numProcs);
//         double end_parallel = MPI_Wtime();
        
//         if (rank == 0) {
//             std::cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << std::endl;
//             std::cout << std::endl;

//             std::ofstream outfile("Parallel_solution.txt");
//             if (!outfile) {
//                 std::cerr << "Error: cannot open file for writing solution.\n";
//                 return -1;
//             }
//             outfile << std::setprecision(15) << std::fixed;
//             for (int i = 0; i < A_csr.n; i++) {
//                 outfile << x[i] << "\n";
//             }
//             outfile.close();
//             std::cout << "Solution x Saved.\n";
//         }
//     } 
//     catch (const std::exception& e) {
//         if (rank == 0)
//             std::cerr << "Error: " << e.what() << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, -1);
//     }
    
//     MPI_Finalize();
//     return 0;
// }




// Hardcoded test case
int main(int argc, char **argv) {
    using namespace std;
    MPI_Init(&argc, &argv);
    
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    try {
        // CSRMatrix L;
        // L.n = 8;
        // L.val = std::vector<double>(14, 1.0); // 
        // L.col_id = {
        //     0,
        //     1,
        //     2,
        //     0, 3,
        //     0, 4,
        //     1, 5,
        //     2, 6,
        //     3, 4, 7
        // };
        // L.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 15};
        // vector<double> b(L.n, 1.0); // b initialized to 1
        // vector<double> x(L.n, 0.0); // x initialized to 0

        CSRMatrix A;
        A.n = 9;
        A.col_id = {
            0,
            1,
            2,
            0, 3,
            0, 4,
            1, 5,
            2, 6,
            3, 4, 7,
            0, 3, 7, 8
        };
        A.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 14, 18};
        A.val = std::vector<double>(18, 1.0);
        vector<double> b(A.n, 1.0); // b initialized to 1
        vector<double> x(A.n, 0.0); // x initialized to 0
        std::vector<std::vector<int>> dependents;
        std::vector<std::vector<int>> dependencies;
        vector<vector<int>> levels = levelScheduling(A, dependents, dependencies);
        if (rank == 0) {
            cout << "Level scheduling result:" << std::endl;
            for (size_t lvl = 0; lvl < levels.size(); lvl++) {
                cout << "Level " << lvl << ": ";
                for (int row : levels[lvl])
                    cout << row << " ";
                cout << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start_parallel = MPI_Wtime();

        parallelTriangularSolve_p2p(
            A, b, x, 
            dependents, dependencies,
            levels, rank, numProcs);
        double end_parallel = MPI_Wtime();
        
        std::vector<double> sendBuf(A.n, 0.0);
        for (int i = 0; i < A.n; ++i) {
            if ((i % numProcs) == rank) {
                sendBuf[i] = x[i];
            }
        }

        std::vector<double> recvBuf(A.n, 0.0);

        // MPI_Reduce，把分散的值加在一起（不会冲突，因为每个行仅在一个rank有非0）
        // root=0，MPI_SUM
        MPI_Reduce(sendBuf.data(), recvBuf.data(), A.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        
        if (rank == 0) {
            cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << std::endl;
            cout << "Parallel solution x:" << std::endl;
            std::cout << "[Rank 0] final solution x:" << std::endl;
            for (int i = 0; i < A.n; ++i) {
                std::cout << recvBuf[i] << " ";
            }
            std::cout << std::endl;
        }
        
    } 
    catch (const exception &e) {
        if (rank == 0)
            cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    MPI_Finalize();
    return 0;
}