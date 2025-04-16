#include <iostream>
#include <vector>
#include <stdexcept>
#include "common.h"

using namespace std;

vector<vector<int>> levelScheduling(const CSRMatrix& A, 
    vector<vector<int>>& dependents,
    vector<vector<int>>& dependencies) 
{
    int n = A.n;
    // for each row j, get all the dependents that relly on j
    vector<int> in_degree(n, 0);
    dependents.assign(n, vector<int>());
    dependencies.assign(n, vector<int>());

    for (int i = 0; i < n; ++i) {  // traverse each row's element 
                                   // idx is the sequence idx of nnz in the row
        for (int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; idx++) {
            int j = A.col_id[idx];
            if (j < i) {
                in_degree[i]++;
                dependents[j].push_back(i);
                dependencies[i].push_back(j);
            }
        }
    }

    vector<vector<int>> levels;   
    vector<int> currentLevel;

    // lev 1：all in-degree 0 nodes
    for (int i = 0; i < n; i++) {
        if (in_degree[i] == 0)
            currentLevel.push_back(i);
    }

    while (!currentLevel.empty()) 
    {
        levels.push_back(currentLevel);
        vector<int> nextLevel;
        // 对当前层的每个节点，将其“删除”（即对依赖它的节点，将入度减1）
        for (int row : currentLevel) {  // row = node
            for (int dependent : dependents[row]) {
                in_degree[dependent]--;
                if (in_degree[dependent] == 0) 
                    nextLevel.push_back(dependent);
            }
        }
        currentLevel = nextLevel;
    }

    // Check for cycles: if any node still has in-degree > 0, there is a cycle
    for (int i = 0; i < n; i++) {
        if (in_degree[i] != 0) {
            throw runtime_error("Cycle detected in dependency graph!");
        }
    }

    return levels;
}


// int main() {
//     CSRMatrix A;
//     A.n = 9;
//     A.col_id = {
//         0,
//         1,
//         2,
//         0, 3,
//         0, 4,
//         1, 5,
//         2, 6,
//         3, 4, 7,
//         0, 3, 7, 8
//     };
//     A.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 14, 19};
//     vector<vector<int>> dependents;
//     vector<vector<int>> dependencies;
//                         // dependencies;
//     try {
//         vector<vector<int>> levels = levelScheduling(A, dependents, dependencies);
//         cout << "Dependency levels:" << endl;
//         for (size_t lvl = 0; lvl < levels.size(); lvl++) {
//             cout << "Level " << lvl << ": ";
//             for (int row : levels[lvl])
//                 cout << row << " ";
//             cout << endl;
//         }
//         cout << "Dependents of each row:" << endl;
//         for (size_t i = 0; i < A.n; i++) {
//             cout << "Row " << i << ": ";
//             for (int dep : dependents[i])
//                 cout << dep << " ";
//             cout << endl;
//         }
//         cout << "Dependencies of each row:" << endl;
//         for (size_t i = 0; i < A.n; i++) {
//             cout << "Row " << i << ": ";
//             for (int in : dependencies[i])
//                 cout << in << " ";
//             cout << endl;
//         }

//     } catch (const exception &e) {
//         cerr << "Error: " << e.what() << endl;
//     }
//     return 0;
// }
