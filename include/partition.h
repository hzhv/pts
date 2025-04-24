#pragma once
#include <vector>

/**
 * @brief Round-robin row partition
 * @param n total number of rows
 * @param P number of processes 
 */
inline std::vector<int> buildRowOwner_mod(int n, int P) {
    std::vector<int> owner(n);
    for (int row = 0; row < n; ++row)
        owner[row] = row % P;
    return owner;
}

/**
 * @brief row partition
 * @param n total number of rows
 * @param P number of processes 
 */
inline std::vector<int> buildRowOwner_block(int n, int P) {
    std::vector<int> owner(n);
    int base = n / P, rem = n % P, start = 0;
    for (int r = 0; r < P; ++r) {
        int len = base + (r < rem ? 1 : 0);
        for (int row = start; row < start + len; ++row) owner[row] = r;
        start += len;
    }
    return owner;
}

/**
 * @brief Block cyclic row partitioning
 * @param n total number of rows
 * @param P number of processes
 * @param B block size, default 256
 */
inline std::vector<int> buildRowOwner_blockCyclic(int n, int P, int B=256) {
    std::vector<int> owner(n);
    for (int row = 0; row < n; ++row)
        owner[row] = (row / B) % P;
    return owner;
}
