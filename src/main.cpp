#include <mpi.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "common.h"

using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    try {
        CSRMatrix L;
        L.n = 8;
        L.val = std::vector<double>(14, 1.0); // 
        L.col_id = {
            0,
            1,
            2,
            0, 3,
            0, 4,
            1, 5,
            2, 6,
            3, 4, 7
        };
        L.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 15};
        vector<double> b(L.n, 1.0); // b initialized to 1
        vector<double> x(L.n, 0.0); // x initialized to 0
        
        vector<vector<int>> levels = levelScheduling(L);
        if (rank == 0) {
            cout << "Level scheduling result:" << endl;
            for (size_t lvl = 0; lvl < levels.size(); lvl++) {
                cout << "Level " << lvl << ": ";
                for (int row : levels[lvl])
                    cout << row << " ";
                cout << endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start_parallel = MPI_Wtime();
        parallelTriangularSolve(L, b, x, levels, rank, numProcs);
        double end_parallel = MPI_Wtime();
        
        if (rank == 0) {
            cout << "Parallel solve time: " << (end_parallel - start_parallel) << " sec" << endl;
            cout << "Parallel solution x:" << endl;
            for (int i = 0; i < L.n; i++) {
                cout << x[i] << " ";
            }
            cout << endl;
        }
        
        // Serial solve for comparison (only on rank 0)
        cout << '\n';
        if (rank == 0) {
            vector<double> x_serial(L.n, 0.0);
            double start_serial = MPI_Wtime();
            serialTriangularSolve(L, b, x_serial);
            double end_serial = MPI_Wtime();
            cout << "Serial solve time: " << (end_serial - start_serial) << " sec" << endl;
            cout << "Serial solution x:" << endl;
            for (int i = 0; i < L.n; i++) {
                cout << x_serial[i] << " ";
            }
            cout << endl;
        }

    } 
    catch (const exception &e) {
        if (rank == 0)
            cerr << "Error: " << e.what() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    MPI_Finalize();
    return 0;
}