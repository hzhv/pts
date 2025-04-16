#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include "../include/common.h"


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

void serialTriangularSolve(const CSRMatrix &L, const std::vector<double>& b, std::vector<double>& x) {
    for (int row = 0; row < L.n; ++row) {
        double sum = 0.0;
        // double diag = 0.0;
        int diag_idx = L.row_ptr[row+1] - 1;
        double diag = L.val[diag_idx];
        for (int idx = L.row_ptr[row]; idx < L.row_ptr[row+1]; ++idx) { // traverse 1 row
            int col = L.col_id[idx];
            sum += L.val[idx] * x[col];
        }
        x[row] = (b[row] - sum) / diag;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <triplet_filename>\n";
        return -1;
    }
    
    std::string CSR_file = argv[1];
    std::ifstream infile(CSR_file);
    std::cout << "Loading A in CSR format..." << std::endl;
    CSRMatrix L = loadCSR(CSR_file);
    std::cout << "CSR Matrix loaded:" << std::endl;
    std::cout << "Number of rows: " << L.n << std::endl;
    std::cout << "Nonzeros: " << L.val.size() << std::endl;
    std::cout << std::endl;
    std::vector<double> b(L.n, 1.0);
    std::vector<double> x(L.n, 0.0);

    std::cout << "Starting serial triangular solve...\n";
    auto start = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 100; ++i)
        serialTriangularSolve(L, b, x);
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    std::cout << "Serial triangular solve time: " << duration << " seconds\n";

    // Save x with higher precision
    std::ofstream outfile("Serial_solution.txt");
    if (!outfile) {
        std::cerr << "Error: cannot open file for writing solution.\n";
        return -1;
    }
    outfile << std::setprecision(15) << std::fixed;
    for (int i = 0; i < L.n; i++) {
        outfile << x[i] << "\n";
    }
    outfile.close();
    std::cout << "Saved.\n";

    return 0;
}