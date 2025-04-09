#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <iomanip>

std::vector<double> loadSolution(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::vector<double> x;
    double value;
    while (infile >> value) {
        x.push_back(value);
    }

    return x;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <serial_solution.txt> <parallel_solution.txt>\n";
        return 1;
    }

    std::string serial_file = argv[1];
    std::string parallel_file = argv[2];

    std::vector<double> x_serial = loadSolution(serial_file);
    std::vector<double> x_parallel = loadSolution(parallel_file);

    if (x_serial.size() != x_parallel.size()) {
        std::cerr << "Error: Solution vectors have different sizes.\n";
        std::cerr << "Serial size: " << x_serial.size() << ", Parallel size: " << x_parallel.size() << std::endl;
        return 1;
    }

    double max_error = 0.0;
    double sum_error = 0.0;
    double tolerance = 1e-12;
    bool passed = true;

    for (size_t i = 0; i < x_serial.size(); ++i) {
        double error = std::abs(x_serial[i] - x_parallel[i]);
        max_error = std::max(max_error, error);
        sum_error += error;
        if (error > tolerance) {
            passed = false;
        }
    }

    double mean_error = sum_error / x_serial.size();

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Max absolute error:  " << max_error << std::endl;
    std::cout << "Mean absolute error: " << mean_error << std::endl;
    std::cout << "Tolerance threshold: " << tolerance << std::endl;
    std::cout << (passed ? "PASS: Solutions match within tolerance." 
                         : "FAIL: Differences exceed tolerance.") << std::endl;

    return 0;
}
