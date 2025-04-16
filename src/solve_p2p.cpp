#include <mpi.h>
#include <vector>
#include <iostream>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "common.h"

using namespace std;


void parallelTriangularSolve_p2p(
    const CSRMatrix& L, 
    const std::vector<double>& b, 
    std::vector<double>& x,
    const std::vector<std::vector<int>>& dependents,
    const std::vector<std::vector<int>>& dependencies,
    const std::vector<std::vector<int>>& levels,
    int rank, 
    int numProcs) 
{
    int n = L.n;

    std::vector<bool> computed(n, false);
    std::vector<bool> received(n, false);
    std::vector<bool> sent(n, false);

    for (size_t lev = 0; lev < levels.size(); ++lev) {

        // Step A: Recv the 
        {
            std::set<int> neededRemotes;
            for (int row : levels[lev]) {
                if ((row % numProcs) == rank && !computed[row]) {
                    for (int dep : dependencies[row]) {
                        if ((dep % numProcs) != rank && !received[dep]) {
                            neededRemotes.insert(dep);
                        }
                    }
                }
            }

            for (int dep : neededRemotes) {
                int sender_rank = dep % numProcs;
                std::cout << "[Rank " << rank << "] ready to recv x[" << dep 
                          << "] from Rank " << sender_rank << " (tag=" << dep << ")" 
                          << std::endl;

                MPI_Recv(&x[dep], 1, MPI_DOUBLE, sender_rank, dep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                received[dep] = true;

                std::cout << "[Rank " << rank << "] received x[" << dep << "]=" << x[dep]
                          << " from Rank " << sender_rank << " (tag=" << dep << ")" 
                          << std::endl;
            }
        }

        // Step B: Calculation
        {
            for (int row : levels[lev]) {
                if ((row % numProcs) == rank && !computed[row]) {
                    double diag = 0.0;
                    double sum  = 0.0;

                    for (int idx = L.row_ptr[row]; idx < L.row_ptr[row + 1]; ++idx) {
                        int col = L.col_id[idx];
                        double val = L.val[idx];
                        if (col == row) {
                            diag = val; 
                        } else {
                            sum += val * x[col];
                        }
                    }
                    x[row] = (b[row] - sum) / diag;
                    computed[row] = true;

                    std::cout << "[Rank " << rank << "] computed x[" << row << "] = " << x[row] 
                              << std::endl;
                }
            }
        }

        // Step C: 将本层新计算好的行发送给更后面层次需要它们的进程
        {
            for (int row : levels[lev]) {
                if ((row % numProcs) == rank && computed[row] && !sent[row]) {
                    std::set<int> target_ranks;
                    for (int dep : dependents[row]) {
                        int dst_rank = dep % numProcs;
                        if (dst_rank != rank) {
                            target_ranks.insert(dst_rank);
                        }
                    }

                    for (int dst : target_ranks) {
                        MPI_Send(&x[row], 1, MPI_DOUBLE, dst, row, MPI_COMM_WORLD);
                        std::cout << "[Rank " << rank << "] sent x[" << row << "]=" << x[row] 
                                  << " to Rank " << dst << " (tag=" << row << ")" << std::endl;
                    }
                    sent[row] = true;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}