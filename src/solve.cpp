#include <mpi.h>
#include <vector>
#include "common.h"

using namespace std;

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

// ------- The following is for serial solving and Testing only -----------------------------
void serialTriangularSolve(const CSRMatrix &L, const vector<double> &b, vector<double> &x) 
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