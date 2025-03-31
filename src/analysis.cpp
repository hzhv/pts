#include <iostream>
#include <vector>
#include <stdexcept>
#include "common.h"

using namespace std;

vector<vector<int>> levelScheduling(const CSRMatrix& A) {
    int n = A.n;

    vector<int> in_degree(n, 0);
    vector<vector<int>> dependents(n);

    for (int i = 0; i < n; i++) {
        for (int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; idx++) {
            int j = A.col_id[idx];
            if (j < i) {
                in_degree[i]++;
                dependents[j].push_back(i);
            }
        }
    }

    vector<vector<int>> levels;  
    vector<int> currentLevel;

    for (int i = 0; i < n; i++) {
        if (in_degree[i] == 0)
            currentLevel.push_back(i);
    }

    while (!currentLevel.empty()) {
        levels.push_back(currentLevel);
        vector<int> nextLevel;
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
