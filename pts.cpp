#include <mpi.h>
#include <vector>
#include <iostream>
#include <cassert>

// CSR 格式存储下三角矩阵
struct CSRMatrix {
    // row_ptr 的长度为 n+1（n 为行数）
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> values;
};

// 消息结构，用于在进程间传递“已计算行”的信息
struct Message {
    int row;       // 已计算行的全局行号
    double value;  // 对应行计算后的解 x[row]
};

// 采用 MPI 的非阻塞通信实现并行下三角求解
// 参数说明：
//   L             —— CSR 格式的下三角矩阵
//   b             —— 右端项向量
//   x             —— 解向量（各进程持有全局解的副本，只有“本地”行会被计算）
//   dependents    —— 对于每一行，给出依赖该行计算结果的所有行号
//   dependencies  —— 对于每一行，给出该行需要等待的上游（依赖）行号
//   levels        —— 分层调度，每一层中行间互不依赖（例如：level[0] 包含所有无依赖的行）
//   rank, numProcs —— 当前进程号及进程总数（行分配采用 row % numProcs 策略）
void parallelTriangularSolve_p2p(
    const CSRMatrix& L, 
    const std::vector<double>& b, 
    std::vector<double>& x,
    const std::vector<std::vector<int>>& dependents,
    const std::vector<std::vector<int>>& dependencies,
    const std::vector<std::vector<int>>& levels,
    int rank, int numProcs)
{
    int n = b.size();
    // 记录每行剩余未满足的依赖数
    std::vector<int> dependency_counter(n, 0);
    // 标记某行是否已经计算完成
    std::vector<bool> processed(n, false);

    // 初始化依赖计数：每行需要等待的依赖数等于 dependencies[row].size()
    for (int i = 0; i < n; i++) {
        dependency_counter[i] = dependencies[i].size();
    }

    // -------------------------------
    // 定义 MPI 消息数据类型 Message
    // -------------------------------
    MPI_Datatype MPI_Message_Type;
    {
        // 为获取地址需要一个 dummy 变量
        Message dummy;
        int block_lengths[2] = {1, 1};
        MPI_Aint displacements[2];
        MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
        MPI_Aint base_address;
        MPI_Get_address(&dummy, &base_address);
        MPI_Get_address(&dummy.row, &displacements[0]);
        MPI_Get_address(&dummy.value, &displacements[1]);
        displacements[0] = displacements[0] - base_address;
        displacements[1] = displacements[1] - base_address;
        MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_Message_Type);
        MPI_Type_commit(&MPI_Message_Type);
    }

    // 用于存储非阻塞发送的请求
    std::vector<MPI_Request> send_requests;

    // 遍历各层，每一层内部的行互不依赖，可并行计算
    for (size_t lvl = 0; lvl < levels.size(); lvl++) {
        const std::vector<int>& current_level = levels[lvl];
        bool level_complete = false;
        while (!level_complete) {

            // 1. 非阻塞轮询检测是否有远程消息到达
            int flag;
            MPI_Status status;
            while (true) {
                MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
                if (!flag)
                    break;
                // 收到消息：获取 Message 信息
                Message msg;
                MPI_Recv(&msg, 1, MPI_Message_Type, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 更新全局解 x 对于 msg.row 的值（各进程均保存一份）
                x[msg.row] = msg.value;
                // 对于依赖于 msg.row 的每个行（依赖由 dependents 给出），若该行归属当前进程，则其依赖计数减 1
                for (int depRow : dependents[msg.row]) {
                    if ((depRow % numProcs) == rank) {
                        dependency_counter[depRow]--;
                    }
                }
            }

            // 2. 处理本层中当前进程负责的、依赖已满足且未处理的行
            bool progress = false;
            for (int row : current_level) {
                // 行分配：采用 row % numProcs 策略
                if ((row % numProcs) != rank)
                    continue;
                if (processed[row])
                    continue;
                if (dependency_counter[row] > 0)
                    continue; // 尚有未满足的依赖

                // 计算该行：
                // 公式： x[row] = (b[row] - sum_{j<row} L(row,j)*x[j]) / L(row,row)
                double sum = 0.0;
                double diag = 0.0;
                for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; idx++) {
                    int col = L.col_idx[idx];
                    double val = L.values[idx];
                    if (col < row)
                        sum += val * x[col];
                    else if (col == row)
                        diag = val;  // 对角元
                }
                x[row] = (b[row] - sum) / diag;
                processed[row] = true;
                progress = true;

                // 3. 将本行的计算结果通知所有依赖于该行的其他行
                for (int depRow : dependents[row]) {
                    if ((depRow % numProcs) == rank) {
                        // 如果依赖行归属当前进程，则直接将其依赖计数减 1
                        dependency_counter[depRow]--;
                    } else {
                        // 否则，构造消息，非阻塞发送至对应进程
                        Message msg;
                        msg.row = row;
                        msg.value = x[row];
                        MPI_Request req;
                        MPI_Isend(&msg, 1, MPI_Message_Type, depRow % numProcs, 0, MPI_COMM_WORLD, &req);
                        send_requests.push_back(req);
                    }
                }
            }

            // 4. 判断本层是否完成：所有本地分配的行均计算完成
            bool allProcessed = true;
            for (int row : current_level) {
                if ((row % numProcs) == rank && !processed[row]) {
                    allProcessed = false;
                    break;
                }
            }
            if (allProcessed)
                level_complete = true;

            // 5. 确保所有非阻塞发送完成（本层完成前同步发送请求）
            if (!send_requests.empty()) {
                MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
                send_requests.clear();
            }
        } // end while 当前层未完成

        // 可选：层与层之间同步，保证所有进程完成当前层后再进入下一层
        MPI_Barrier(MPI_COMM_WORLD);
    } // end for 每一层

    // 释放自定义 MPI 数据类型
    MPI_Type_free(&MPI_Message_Type);
}


// ========================
// 主函数：构造示例系统并调用求解函数
// ========================
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // 构造一个简单的 4×4 下三角矩阵
    // 示例矩阵（下三角部分非零）：
    //   L = [ 2         0    0   0 ]
    //       [ 3         4    0   0 ]
    //       [ 0         1    5   0 ]
    //       [ 2         0    3   6 ]
    // 对应 CSR 格式：
    //   row_ptr:  {0, 1, 3, 5, 8}
    //   col_idx:  {0, 0, 1, 1, 2, 0, 2, 3}
    //   values:   {2.0, 3.0, 4.0, 1.0, 5.0, 2.0, 3.0, 6.0}
    CSRMatrix L;
    L.row_ptr = {0, 1, 3, 5, 8};
    L.col_idx = {0, 0, 1, 1, 2, 0, 2, 3};
    L.values  = {2.0, 3.0, 4.0, 1.0, 5.0, 2.0, 3.0, 6.0};

    int n = 4;
    // 右端项向量 b
    std::vector<double> b = {2.0, 7.0, 8.0, 18.0};
    // 初始化解向量 x（所有元素置 0）
    std::vector<double> x(n, 0.0);

    // 构造 dependencies（每行依赖于哪些上游行）
    // 行 0：无依赖
    // 行 1：依赖行 0（由于 L(1,0) ≠ 0）
    // 行 2：依赖行 1（由于 L(2,1) ≠ 0）
    // 行 3：依赖行 0 和 2（由于 L(3,0) 与 L(3,2) ≠ 0）
    std::vector<std::vector<int>> dependencies(n);
    dependencies[0] = {};
    dependencies[1] = {0};
    dependencies[2] = {1};
    dependencies[3] = {0, 2};

    // 构造 dependents（依赖于该行结果的其他行）
    // 行 0：依赖行有 1 和 3
    // 行 1：依赖行有 2
    // 行 2：依赖行有 3
    // 行 3：无依赖
    std::vector<std::vector<int>> dependents(n);
    dependents[0] = {1, 3};
    dependents[1] = {2};
    dependents[2] = {3};
    dependents[3] = {};

    // 根据依赖关系，构造层次调度（level scheduling）
    // 此处采用简单层次：先计算无依赖的行，再依次计算其依赖行
    // level 0: {0}
    // level 1: {1}
    // level 2: {2}
    // level 3: {3}
    std::vector<std::vector<int>> levels;
    levels.push_back({0});
    levels.push_back({1});
    levels.push_back({2});
    levels.push_back({3});

    // 调用并行三角求解函数
    parallelTriangularSolve_p2p(L, b, x, dependents, dependencies, levels, rank, numProcs);

    // 为了便于展示，每个进程都持有一份解向量，但可能只有部分行为本地计算结果
    // 使用 MPI_Allreduce 聚合所有进程计算得到的 x（取最大值即可，由于各进程对同一行计算结果一致）
    std::vector<double> x_global(n, 0.0);
    MPI_Allreduce(x.data(), x_global.data(), n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // 由 rank 0 输出最终解向量
    if (rank == 0) {
        std::cout << "求解结果 x:" << std::endl;
        for (int i = 0; i < n; i++) {
            std::cout << "x[" << i << "] = " << x_global[i] << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}