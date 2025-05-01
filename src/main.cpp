#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "common.h"
#include "partition.h"

// ------------------------------------------------------------------------------
template<typename T>
void bcast_vector(std::vector<T>& vec, int root, MPI_Comm comm)
{
    int sz = static_cast<int>(vec.size());
    MPI_Bcast(&sz, 1, MPI_INT, root, comm);
    if (vec.empty()) 
        vec.resize(sz);
    MPI_Datatype dtype = (std::is_same<T,int>::value ? MPI_INT : MPI_DOUBLE);
    MPI_Bcast(vec.data(), sz, dtype, root, comm);
}
// ------------------------------------------------------------------------------

int main(int argc, char* argv[]) 
{
    MPI_Init(&argc, &argv);
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (argc < 2 && rank == 0) 
    {
        std::cerr << "Usage: " << argv[0] << " <csr_filename>\n";
    }
    if (argc < 2) 
    { 
        MPI_Finalize(); 
        return -1; 
    }

    // Get and Bcast A_csr
    CSRMatrix A_csr;
    if (rank == 0) 
    {
        std::string CSR_file = argv[1];
        std::ifstream infile(CSR_file);
        if (!infile) {
            std::cerr << "Cannot open file: " << CSR_file << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        std::cout << "Loading A in CSR format..." << std::endl;
        A_csr = loadCSR(CSR_file);
        std::cout << "[rank0]  Matrix loaded: n = "
                  << A_csr.n << ", nnz = " << A_csr.val.size() << '\n';
    }

    MPI_Bcast(&A_csr.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    bcast_vector(A_csr.row_ptr, 0, MPI_COMM_WORLD);
    bcast_vector(A_csr.col_id,  0, MPI_COMM_WORLD);
    bcast_vector(A_csr.val,     0, MPI_COMM_WORLD);

    // Get and Bcast x and b
    std::vector<double> b;
    std::vector<double> x(A_csr.n, 0.0);
    if (rank == 0) 
    {
        b.assign(A_csr.n, 1.0);
        x.assign(A_csr.n, 0.0);
    }
    bcast_vector(b, 0, MPI_COMM_WORLD);
    bcast_vector(x, 0, MPI_COMM_WORLD);

    if (rank == 0) 
    {
        std::cout << "Broadcast of CSRMatrix completed. Matrix n = " << A_csr.n << std::endl;
        std::cout << std::endl;
    }

    std::vector<int> row_owner = buildRowOwner_mod(A_csr.n, numProcs);
    // ====================== "Fastest" Test start ======================

    std::vector<int> level_ptr, level_rows,
                     level_counts_flat, level_displs_flat,
                     send_rows_ptr, send_rows,
                     recv_perm_ptr,  recv_perm; // lD Array

    buildSchedules
    (
        A_csr, row_owner, rank, numProcs,
        level_ptr, level_rows,
        level_counts_flat, level_displs_flat,
        send_rows_ptr, send_rows,
        recv_perm_ptr, recv_perm
    );

    double comm_time = 0.0;
    double start_parallel = MPI_Wtime();
    try 
    {   
        for (int i = 0; i < 1000; ++i)  
            parallelTriangularSolve_fast
            (
                A_csr, level_ptr, level_rows,
                level_counts_flat, level_displs_flat,
                send_rows_ptr, send_rows,
                recv_perm_ptr,  recv_perm,
                b, x, 
                comm_time,
                rank, numProcs
            );
    }
    // ====================== "Fastest" Test End ======================
    // std::vector<int> level_ptr, level_rows, dep_ptr, dep_rows;
    // if (rank == 0) 
    // {   
    //     levelScheduling_plain(A_csr, level_ptr, level_rows, dep_ptr, dep_rows);
    // }
    // bcast_vector(level_ptr,  0, MPI_COMM_WORLD);
    // bcast_vector(level_rows, 0, MPI_COMM_WORLD);    
    // if (rank == 0) 
    // {
    //     std::cout << "Broadcast of levels and dependents completed.\n" << 
    //     "Levels 1D array length = " << level_ptr.size()-1 << std::endl; // = # of lvls
    //     // "\nDependents 1D array size = " << dep_ptr.size()-1 << std::endl; // = # of rows
    // }

    // double start_parallel = MPI_Wtime();
    // try 
    // {   
    //     for (int i = 0; i < 1000; ++i)
    //         parallelTriangularSolve_block
    //         (
    //             A_csr, 
    //             x, b, 
    //             level_ptr, 
    //             level_rows,
    //             row_owner,
    //             dep_ptr, dep_rows,
    //             rank, numProcs
    //         );
    // } 
    // ====================== Block Test End ======================
    catch (const std::exception& e) 
    {
        if (rank == 0)
        std::cerr << "Error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    double end_parallel = MPI_Wtime();
    
    if (rank == 0) 
    {
        // std::cout << "Total communication time: " 
        //           << comm_time/1000.0 << " sec" << "\n";
        std::cout << "Parallel solve time: " 
            << (end_parallel - start_parallel) / 1000.0 << " sec" << "\n";

        std::ofstream outfile("Parallel_solution.txt");
        if (!outfile) {
            std::cerr << "Error: cannot open file for writing solution.\n";
            return -1;
        }
        outfile << std::setprecision(16) << std::fixed;
        for (int i = 0; i < A_csr.n; i++) {
            outfile << x[i] << "\n";
        }
        outfile.close();
        std::cout << "Solution x Saved.\n";
    }

    MPI_Finalize();
    return 0;
}

// // pts_p2p Test code:
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
//         MPI_Barrier(MPI_COMM_WORLD);
//         double start_parallel = MPI_Wtime();
//         parallelTriangularSolve_p2p(
//             A_csr, b, x, 
//             dependents, dependencies,
//             levels, rank, numProcs);
//         double end_parallel = MPI_Wtime();
        
//         std::vector<double> sendBuf(A_csr.n, 0.0);
//         for (int i = 0; i < A_csr.n; ++i) {
//             if ((i % numProcs) == rank) {
//                 sendBuf[i] = x[i];
//             }
//         }
//         std::vector<double> recvBuf;
//         if (rank == 0 ) recvBuf.resize(A_csr.n, 0.0);

//         // MPI_Reduce，把分散的值加在一起（不会冲突，因为每个行仅在一个rank有非0）
//         MPI_Reduce(sendBuf.data(), 
//                    rank == 0 ? recvBuf.data() : nullptr,
//                    A_csr.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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




// // Hardcoded test case
// int main(int argc, char **argv) {
//     using namespace std;
//     MPI_Init(&argc, &argv);
    
//     int rank, numProcs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

//     try {
//         // CSRMatrix L;
//         // L.n = 8;
//         // L.val = std::vector<double>(14, 1.0); // 
//         // L.col_id = {
//         //     0,
//         //     1,
//         //     2,
//         //     0, 3,
//         //     0, 4,
//         //     1, 5,
//         //     2, 6,
//         //     3, 4, 7
//         // };
//         // L.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 15};
//         // vector<double> b(L.n, 1.0); // b initialized to 1
//         // vector<double> x(L.n, 0.0); // x initialized to 0

        // CSRMatrix A;
        // A.n = 9;
        // A.col_id = {
        //     0,
        //     1,
        //     2,
        //     0, 3,
        //     0, 4,
        //     1, 5,
        //     2, 6,
        //     3, 4, 7,
        //     0, 3, 7, 8
        // };
        // A.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 14, 18};
        // A.val = std::vector<double>(18, 1.0);
//         vector<double> b(A.n, 1.0); // b initialized to 1
//         vector<double> x(A.n, 0.0); // x initialized to 0
//         std::vector<std::vector<int>> dependents;
//         std::vector<std::vector<int>> dependencies;
//         vector<vector<int>> levels = levelScheduling(A, dependents, dependencies);
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
//         parallelTriangularSolve_p2p(
//             A, b, x, 
//             dependents, dependencies,
//             levels, rank, numProcs);
//         double end_parallel = MPI_Wtime();
        
//         std::vector<double> sendBuf(A.n, 0.0);
//         for (int i = 0; i < A.n; ++i) {
//             if ((i % numProcs) == rank) {
//                 sendBuf[i] = x[i];
//             }
//         }
//         std::vector<double> recvBuf;
//         if (rank == 0 ) recvBuf.resize(A.n, 0.0);

//         // MPI_Reduce，把分散的值加在一起（不会冲突，因为每个行仅在一个rank有非0）
//         MPI_Reduce(sendBuf.data(), 
//                    rank == 0 ? recvBuf.data() : nullptr,
//                    A.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//         if (rank == 0) {
//             cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << std::endl;
//             cout << "Parallel solution x:" << std::endl;
//             std::cout << "[Rank 0] final solution x:" << std::endl;
//             for (int i = 0; i < A.n; ++i) {
//                 std::cout << recvBuf[i] << " ";
//             }
//             std::cout << std::endl;
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