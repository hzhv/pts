#include <iostream>
#include <vector>
#include <stdexcept>
#include "common.h"

using namespace std;

// 使用level scheduling对稀疏下三角矩阵构造依赖分层
// 返回一个二维数组levels，其中levels[k]包含第k层(从0开始)可以并行计算的行索引
vector<vector<int>> levelScheduling(const CSRMatrix& A) {
    int n = A.n;
    // In-degree: for row i, gather all j < i and A(i,j)!=0)
    vector<int> in_degree(n, 0);
    // for each row j, get all the dependents that relly on j
    vector<vector<int>> dependents(n);

    for (int i = 0; i < n; i++) {
        for (int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; idx++) {
            int j = A.col_id[idx];
            // 对于下三角矩阵，只考虑 j < i 的情况
            if (j < i) {
                in_degree[i]++;
                dependents[j].push_back(i);
            }
        }
    }

    vector<vector<int>> levels;   
    vector<int> currentLevel;

    // 第一层：所有入度为0的节点
    for (int i = 0; i < n; i++) {
        if (in_degree[i] == 0)
            currentLevel.push_back(i);
    }

    // 反复处理，直到所有节点都被分层
    while (!currentLevel.empty()) {
        levels.push_back(currentLevel);
        vector<int> nextLevel;
        // 对当前层的每个节点，将其“删除”（即对依赖它的节点，将入度减1）
        for (int node : currentLevel) {
            for (int dependent : dependents[node]) {
                in_degree[dependent]--;
                if (in_degree[dependent] == 0) {
                    nextLevel.push_back(dependent);
                }
            }
        }
        currentLevel = nextLevel;
    }

    // 检查是否所有节点都处理完（无环图时应当如此）
    for (int i = 0; i < n; i++) {
        if (in_degree[i] != 0) {
            throw runtime_error("Cycle detected in dependency graph!");
        }
    }

    return levels;
}


// int main() {
//     CSRMatrix A;
//     A.n = 8;
//     A.col_id = {
//         0,
//         1,
//         2,
//         0, 3,
//         0, 4,
//         1, 5,
//         2, 6,
//         3, 4, 7
//     };
//     A.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 15};

//     try {
//         vector<vector<int>> levels = levelScheduling(A);
//         cout << "Dependency levels:" << endl;
//         for (size_t lvl = 0; lvl < levels.size(); lvl++) {
//             cout << "Level " << lvl << ": ";
//             for (int row : levels[lvl])
//                 cout << row << " ";
//             cout << endl;
//         }
//     } catch (const exception &e) {
//         cerr << "Error: " << e.what() << endl;
//     }
//     return 0;
// }
