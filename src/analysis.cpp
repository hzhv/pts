#include <iostream>
#include <vector>
#include <stdexcept>
#include "common.h"

using namespace std;

vector<vector<int>> levelScheduling
/* @p
*/
(
    const CSRMatrix& A, 
    vector<vector<int>>& dependents,
    vector<vector<int>>& dependencies
) 
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
        for (int row : currentLevel) {  // row = node
            for (int dependent : dependents[row]) {
                in_degree[dependent]--;
                if (in_degree[dependent] == 0) 
                    nextLevel.push_back(dependent);
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

void levelScheduling_plain
(
    const CSRMatrix& A,
    std::vector<int>& level_ptr,
    std::vector<int>& level_rows,
    std::vector<int>& dep_ptr,
    std::vector<int>& dep_rows
) 
{
    int n = A.n;
    std::vector<int> in_deg(n, 0);
    std::vector<int> dep_count(n, 0);
    for (int i = 0; i < n; ++i) 
    {
        for (int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; ++idx) 
        {
            int j = A.col_id[idx];
            if (j < i) 
            {
                in_deg[i]++;      
                dep_count[j]++;   
            }
        }
    }

    dep_ptr.resize(n+1);
    dep_ptr[0] = 0;
    for (int i = 0; i < n; ++i) 
    {
        dep_ptr[i+1] = dep_ptr[i] + dep_count[i];
    }
    int m_dep = dep_ptr[n];
    dep_rows.resize(m_dep);

    std::vector<int> cursor = dep_ptr;
    for (int i = 0; i < n; ++i) 
    {
        for (int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; ++idx) 
        {
            int j = A.col_id[idx];
            if (j < i) 
            {
                dep_rows[cursor[j]++] = i;
            }
        }
    }

    level_ptr.clear();
    level_rows.clear();
    level_ptr.push_back(0);

    std::vector<int> curr;
    curr.reserve(n);
    for (int i = 0; i < n; ++i) // O(n) 
    {
        if (in_deg[i] == 0)
            curr.push_back(i);
    }

    // Kahn
    while (!curr.empty()) // O()
    {
        for (int r : curr)
            level_rows.push_back(r);
        level_ptr.push_back((int)level_rows.size());

        std::vector<int> next;
        next.reserve(curr.size()); // next lvl size ussually smaller than curr lvl size
        for (int r : curr) {
            for (int p = dep_ptr[r]; p < dep_ptr[r+1]; ++p) {
                int child = dep_rows[p];
                if (--in_deg[child] == 0) {
                    next.push_back(child);
                }
            }
        }
        curr.swap(next);
    }

    for (int i = 0; i < n; ++i) {
        if (in_deg[i] != 0) {
            throw std::runtime_error("Cycle detected in dependency graph!");
        }
    }
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
//     A.row_ptr = {0, 1, 2, 3, 5, 7, 9, 11, 14, 18};
//     vector<vector<int>> dependents;
//     vector<vector<int>> dependencies;
//                         // dependencies;
    
//     try {
        // vector<vector<int>> levels = levelScheduling(A, dependents, dependencies);
        // cout << "Dependency levels:" << endl;
        // for (size_t lvl = 0; lvl < level.size(); lvl++) {
        //     cout << "Level " << lvl << ": ";
        //     for (int row : levels[lvl])
        //         cout << row << " ";
        //     cout << endl;
        // }
        // cout << "Dependents of each row:" << endl;
        // for (size_t i = 0; i < A.n; i++) {
        //     cout << "Row " << i << ": ";
        //     for (int dep : dependents[i])
        //         cout << dep << " ";
        //     cout << endl;
        // }
        // cout << "Dependencies of each row:" << endl;
        // for (size_t i = 0; i < A.n; i++) {
        //     cout << "Row " << i << ": ";
        //     for (int in : dependencies[i])
        //         cout << in << " ";
        //     cout << endl;
        // }
        
//         std::vector<int> level_ptr, level_rows, dep_ptr, dep_rows;
//         levelSchedulingCSR(A, level_ptr, level_rows, dep_ptr, dep_rows);
//         std::cout << "===== Levels =====\n";
//         for (size_t k = 0; k + 1 < level_ptr.size(); ++k) {
//             std::cout << "Level " << k << ": ";
//             for (int p = level_ptr[k]; p < level_ptr[k + 1]; ++p)
//                 std::cout << level_rows[p] << ' ';
//             std::cout << '\n';
//         }

//         std::cout << "\n===== Dependents (row -> list) =====\n";
//         for (int i = 0; i < A.n; ++i) {
//             std::cout << i << " : ";
//             for (int p = dep_ptr[i]; p < dep_ptr[i + 1]; ++p)
//                 std::cout << dep_rows[p] << ' ';
//             std::cout << '\n';
//         }


//     } catch (const exception &e) {
//         cerr << "Error: " << e.what() << endl;
//     }
//     return 0;
// }



// ==============================================
// CSRMatrix A; 
// A.n = 9; 
// A.col_id = 
// { 
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

// row_ptr_T = {0,4,6,8,11,13,14,15,17,18};
// col_id_T  = {
//     0,3,4,8, // row 0
//     1,5,     // row 1
//     2,6,     // row 2
//     3,7,8,   // row 3
//     4,7,     // row 4
//     5,       // row 5
//     6,       // row 6
//     7,8,     // row 7
//     8        // row 8
// };
