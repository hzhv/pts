#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// 假设 n 为全局矩阵大小，L 为 n*n 下三角矩阵，b 为 n 维向量，x 为解向量
// 假设数据已经初始化，并且每个进程仅保存自己负责的行

double** allocate_matrix(int rows, int cols) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*) malloc(cols * sizeof(double));
    }
    return matrix;
}

double* allocate_vector(int size) {
    return (double*) malloc(size * sizeof(double));
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 100; // n x n 
    int rowsPerProc = n / size; 
    int remainder = n % size;
    
    /* Get the rows of current process deal with: [start, end)
     * rank must smaller than np
     * e.g. n = 101; np = 5; rank = 3
     * rowsPerProc = 100 / 5 = 20
     * remainder = 101 % 5 = 1
     * start = 3 * 20 + (3 < 1 ? 3 : 1) = 60 + 1 = 61
     * local_rows = rowsPerProc + (3 < 1 ? 1 : 0) = 20 + 0 = 20
     * end = start + local_rows = 61 + 20 = 81
     */
    int start = rank * rowsPerProc + (rank < remainder ? rank : remainder);
    int local_rows = rowsPerProc + (rank < remainder ? 1 : 0);
    int end = start + local_rows;
    
    // 分配本地存储：L_local[local_rows][n], b_local[local_rows], x[全局维度]
    double** L_local = allocate_matrix(local_rows, n);
    double* b_local = allocate_vector(local_rows);
    double* x = (double*) malloc(n * sizeof(double));  // 假设每个进程保存全局解（或可用全局缓冲区）
    
    // 初始化 x 数组为 0 或者未定义
    for (int i = 0; i < n; i++) x[i] = 0.0;

    // 主循环：按全局行号迭代
    for (int i = 0; i < n; i++) {
        // 如果当前行 i 不在本进程负责区，则等待接收
        if (i < start || i >= end) {
            // 如果后续有进程需要该 x[i]，本进程可能需要接收它
            MPI_Bcast(&x[i], 1, MPI_DOUBLE, owner_of(i), MPI_COMM_WORLD);
        } 
        else {
            // i 属于本进程：计算对应局部行 index = i - start
            int local_i = i - start;
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L_local[local_i][j] * x[j];
            }
            x[i] = (b_local[local_i] - sum) / L_local[local_i][i];
            // 计算完成后，将 x[i] 广播给所有进程
            MPI_Bcast(&x[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        }
    }

    // 现在每个进程拥有全局解向量 x 中自己负责的部分
    // 可选：收集所有结果到根进程
    // MPI_Gather(...)

    MPI_Finalize();
    return 0;
}
