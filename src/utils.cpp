#include <algorithm>
#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "common.h"


CSRMatrix loadCSR(std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open CSR file: " + filename);
    }
    auto start = std::chrono::high_resolution_clock::now();

    CSRMatrix csr;
    int row_ptr_size, col_id_size, val_size;

    infile >> csr.n;
    infile >> row_ptr_size;
    csr.row_ptr.resize(row_ptr_size);
    for (int i = 0; i < row_ptr_size; ++i) {
        infile >> csr.row_ptr[i];
    }

    infile >> col_id_size;
    csr.col_id.resize(col_id_size);
    for (int i = 0; i < col_id_size; ++i) {
        infile >> csr.col_id[i];
    }

    infile >> val_size;
    csr.val.resize(val_size);
    for (int i = 0; i < val_size; ++i) {
        infile >> csr.val[i];
    }
    infile.close();

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    std::cout << "CSR matrix loaded from \"" << filename << "\" in " 
              << duration << " seconds." << std::endl;
              
    return csr;
}

// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <triplet_filename>\n";
//         return -1;
//     }

//     std::string filename = argv[1];
//     CSRMatrix csr = loadCSR(filename);
//     std::cout << csr.n << std::endl;
//     std::cout << csr.row_ptr.size() << std::endl;
//     std::cout << csr.col_id.size() << std::endl;
//     std::cout << csr.val.size() << std::endl;
//     return 0;
    
// }