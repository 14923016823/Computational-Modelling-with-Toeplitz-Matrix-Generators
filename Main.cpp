#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <chrono>

#include "Vectord.h"
#include "SparseToeplitz.h"
#include "BlockToeplitz.h"

using namespace std::chrono;

// Benchmark function
template<typename Func>
double benchmark(Func f, int iterations = 10) {
    auto start = high_resolution_clock::now();
    for(int i = 0; i < iterations; i++) {
        f();
    }
    auto end = high_resolution_clock::now();
    return duration_cast<microseconds>(end - start).count() / (1000.0 * iterations);
}

int main(int argc, char* argv[])
{
    std::cout << "SparseToeplitz Storage Example:\n";
    std::cout << "For a 4x4 periodic Laplacian matrix:\n";
    std::cout << "Matrix = [2, -1,  0, -1;\n";
    std::cout << "         -1,  2, -1,  0;\n";
    std::cout << "          0, -1,  2, -1;\n";
    std::cout << "         -1,  0, -1,  2]\n\n";
    
    std::cout << "SparseToeplitz stores ONLY:\n";
    std::cout << "Diags = [0, 1, 3] (diagonal offsets)\n";
    std::cout << "Vals  = [2, -1, -1] (values for those diagonals)\n";
    std::cout << "Total storage: 3 integers + 3 doubles = O(1) regardless of matrix size!\n\n";
    
    // Test with small matrix
    int n = 4;
    std::vector<int> diags = {0, 1, 3};
    std::vector<double> vals = {2.0, -1.0, -1.0};
    
    SparseToeplitz<double> toeplitz(n, n, 3, diags.data(), vals.data());
    
    std::cout << "Actual matrix representation:\n";
    toeplitz.print();
    std::cout << "\n";
    
    std::cout << "Benchmarking matvec implementations...\n\n";
    
    // Test sizes
    std::vector<int> sizes = {32, 64, 128, 256, 512};
    
    for(int n : sizes) {
        std::cout << "Matrix size: " << n << "x" << n << "\n";
        
        // Create a circulant Toeplitz matrix (like for periodic BC Laplacian)
        // Only stores non-zero diagonals: main (2), +1 (-1), -1 (-1 which wraps to n-1)
        // Storage: Diags = [0, 1, n-1], Vals = [2.0, -1.0, -1.0] - NO zeros stored!
        std::vector<int> diags = {0, 1, n-1};
        std::vector<double> vals = {2.0, -1.0, -1.0};
        
        SparseToeplitz<double> toeplitz(n, n, 3, diags.data(), vals.data());
        
        // Random input vector
        Vectord input(n);
        for(int i = 0; i < n; i++) {
            input[i] = (double)rand() / RAND_MAX;
        }
        
        Vectord output_regular(n);
        Vectord output_fft(n);
        
        // Benchmark regular matvec
        double time_regular = benchmark([&]() {
            toeplitz.matvec(input, output_regular);
        });
        
        // Benchmark FFT matvec
        double time_fft = benchmark([&]() {
            toeplitz.matvec_fft(input, output_fft);
        });
        
        // Check accuracy
        double max_diff = 0.0;
        for(int i = 0; i < n; i++) {
            max_diff = std::max(max_diff, std::abs(output_regular[i] - output_fft[i]));
        }
        
        std::cout << "  Regular matvec: " << time_regular << " ms\n";
        std::cout << "  FFT matvec:     " << time_fft << " ms\n";
        std::cout << "  Speedup:        " << time_regular / time_fft << "x\n";
        std::cout << "  Max difference: " << max_diff << "\n\n";
    }
    
    // Test BlockToeplitz
    std::cout << "BlockToeplitz benchmark (4x4 blocks of 32x32):\n";
    int block_size = 32;
    int num_blocks = 4;
    int total_size = block_size * num_blocks;
    
    // Create block Toeplitz with 3 diagonals
    BlockToeplitz<double> block_toeplitz(total_size, total_size, 3);
    
    // Set up the diagonals with SparseToeplitz blocks
    for(int d = 0; d < 3; d++) {
        std::vector<int> diags = {0, 1, block_size-1};
        std::vector<double> vals = {2.0, -1.0, -1.0};
        SparseToeplitz<double>* block = new SparseToeplitz<double>(block_size, block_size, 3, diags.data(), vals.data());
        block_toeplitz.set_diag(d, d-1, block); // diagonals at -1, 0, 1
    }
    
    Vectord block_input(total_size);
    Vectord block_output_regular(total_size);
    Vectord block_output_fft(total_size);
    
    for(int i = 0; i < total_size; i++) {
        block_input[i] = (double)rand() / RAND_MAX;
    }
    
    // Benchmark block matvec
    double time_block_regular = benchmark([&]() {
        block_toeplitz.regular_matvec(block_input, block_output_regular);
    });
    
    double time_block_fft = benchmark([&]() {
        block_toeplitz.matvec(block_input, block_output_fft);
    });
    
    // Check accuracy
    double max_diff_block = 0.0;
    for(int i = 0; i < total_size; i++) {
        max_diff_block = std::max(max_diff_block, std::abs(block_output_regular[i] - block_output_fft[i]));
    }
    
    std::cout << "  Block regular matvec: " << time_block_regular << " ms\n";
    std::cout << "  Block FFT matvec:     " << time_block_fft << " ms\n";
    std::cout << "  Speedup:              " << time_block_regular / time_block_fft << "x\n";
    std::cout << "  Max difference:       " << max_diff_block << "\n";
    
    std::cout << "\nBenchmark completed.\n";
    return 0;
}
