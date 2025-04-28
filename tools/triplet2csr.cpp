#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "../include/common.h"

CSRMatrix loadCSRFromTripletFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    std::vector<Triplet> triplets;
    std::string line;
    
    while (std::getline(infile, line)) { // line by line
        if (line.empty()) continue;
        std::istringstream iss(line);
        int row, col;
        double value;

        if (!(iss >> row >> col >> value)) {
            continue;
        }
        // Due to MATLAB's row is 1-based，change to 0-based
        triplets.push_back({row - 1, col, value});
    }
    infile.close();

    int n = 0;
    for (const auto& t : triplets) {
        if (t.row >= n) {
            n = t.row + 1;
        }
    }

    std::sort(triplets.begin(), triplets.end(), [](const Triplet& a, const Triplet& b) 
    {
        return (a.row == b.row) ? (a.col < b.col) : (a.row < b.row);
    });

    int nnz = triplets.size();
    std::vector<int> row_ptr(n + 1, 0);
    std::vector<int> col_id(nnz);
    std::vector<double> values(nnz);

    for (const auto& t : triplets) {
        row_ptr[t.row + 1]++;
    }

    for (int i = 0; i < n; ++i) {
        row_ptr[i + 1] += row_ptr[i];
    }

    std::vector<int> current_row_ptr = row_ptr;

    for (const auto& t : triplets) {
        int pos = current_row_ptr[t.row]++;
        col_id[pos] = t.col;
        values[pos] = t.val;
    }

    CSRMatrix csr;
    csr.n = n;
    csr.row_ptr = std::move(row_ptr);
    csr.col_id = std::move(col_id);
    csr.val = std::move(values);
    return csr;
}

void saveCSRToFile(const CSRMatrix& csr, const std::string& filename) {
    /* Output format:
     * n
     * row_ptr.size()
     * row_ptr
     * col_id.size()
     * col_id
     * val.size()
     * val
     */
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    out << csr.n << "\n"; 
    out << csr.row_ptr.size() << "\n";
    for (int val : csr.row_ptr) {
        out << val << " ";
    }
    out << "\n";

    out << csr.col_id.size() << "\n";
    for (int col : csr.col_id) {
        out << col << " ";
    }
    out << "\n";

    out << csr.val.size() << "\n";
    for (double v : csr.val) {
        out << std::setprecision(16) << v << " ";
    }
    out << "\n";

    out.close();
    std::cout << "CSR saved to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <triplet_file>\n";
        return -1;
    }

    std::string triplet_filename = argv[1];
    std::string output_csr_file = triplet_filename + "_csr.txt";  

    std::cout << "Loading triplet matrix..." << std::endl;
    CSRMatrix A = loadCSRFromTripletFile(triplet_filename);
    std::cout << "Saving CSR..." << std::endl;
    saveCSRToFile(A, output_csr_file);

    return 0;
}

// int main() {
//     std::ifstream fin("L_triplet.txt");
//     if (!fin.is_open()) {
//         std::cerr << "Cannot open file.\n";
//         return 1;
//     } else {
//         std::cout << "Loading... \n";
//     }

//     std::vector<Triplet> triplets;
//     std::string line;

//     while (std::getline(fin, line)) {
//         std::istringstream iss(line);
//         Triplet t;
//         if (iss >> t.row >> t.col >> t.val) {
//             triplets.push_back(t);
//         } else {
//             std::cerr << "Skipping invalid line: " << line << "\n";
//         }
//     }

//     fin.close();

//     std::cout << "Loaded " << triplets.size() << " entries.\n";
//     std::cout.precision(16);
//     for (int i = 0; i < std::min(5, (int)triplets.size()); i++) {
//         std::cout << triplets[i].row << " "
//                   << triplets[i].col << " "
//                   << triplets[i].val << "\n";
//     }

//     return 0;
// }
