#include <iostream>
#include <vector>
#include <exception>
#include "cmmcore/memoryview.h"
// Assume the MemoryView and Slice classes are defined as above
using namespace cmmcore;
#include <iostream>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <typeinfo>
#include <cstring>
#include <memory>
#include <cstdio>

// Assume the enhanced MemoryView class is defined as above
int test2() {
    try {
        // 1. Create a 3D array with shape (3, 4, 5)
        std::vector<int> data(3 * 4 * 5);
        // Initialize data with sequential values for demonstration
        std::iota(data.begin(), data.end(), 0);

        MemoryView<int> mem1(data.data(), {3, 4, 5});
        std::cout << "Original MemoryView:\n" << mem1 << "\n\n";

        // 2. Perform slicing: mem1.slice({Slice(), Slice(1, 3), Slice()})
        // Equivalent to mem1[:, 1:3, :]
        std::vector<Slice> slices = { Slice(), Slice(1, 3), Slice() };
        MemoryView<int> mem2 = mem1.slice(slices);
        std::cout << "After slicing with {{}, {1, 3}, {}}:\n" << mem2 << "\n\n";

        // 3. Reshape mem2 to (3, 2, 5)
        mem2.reshape({3, 2, 5});
        std::cout << "After reshaping mem2 to (3, 2, 5):\n" << mem2 << "\n\n";

        // 4. Access an element: mem2( {2, 1, 4} )
        int value = mem2({2, 1, 4});
        std::cout << "Element at (2, 1, 4): " << value << "\n\n";

        // 5. Serialize mem2
        std::vector<char> serialized_data = mem2.serialize();
        std::cout << "Serialized mem2 into buffer of size " << serialized_data.size() << " bytes.\n\n";

        // 6. Deserialize to create mem3
        MemoryView<int> mem3 = MemoryView<int>::deserialize(serialized_data);
        std::cout << "Deserialized MemoryView (mem3):\n" << mem3 << "\n\n";

        // 7. Write mem3 to a file
        FILE* outfile = fopen("memoryview.dat", "wb");
        if (!outfile) {
            throw std::runtime_error("Failed to open file for writing");
        }
        mem3.write_to_file(outfile);
        fclose(outfile);
        std::cout << "Serialized mem3 written to 'memoryview.dat'.\n\n";

        // 8. Read from the file to create mem4
        FILE* infile = fopen("memoryview.dat", "rb");
        if (!infile) {
            throw std::runtime_error("Failed to open file for reading");
        }
        MemoryView<int> mem4 = MemoryView<int>::read_from_file(infile);
        fclose(infile);
        std::cout << "Deserialized MemoryView from 'memoryview.dat' (mem4):\n" << mem4 << "\n\n";

        // 9. Attempt to deserialize with type mismatch
        try {
            // Suppose we try to deserialize as double instead of int
            std::vector<char> serialized_double = mem2.serialize(); // Still int type
            MemoryView<double> mem_invalid = MemoryView<double>::deserialize(serialized_double);
        } catch (const std::exception& e) {
            std::cout << "Caught exception for type mismatch during deserialization: " << e.what() << "\n\n";
        }

        // 10. Attempt to reshape to incompatible shape
        try {
            mem2.reshape({2, 2, 5}); // Original total elements: 30
            // New shape total elements: 20 (invalid)
        } catch (const std::exception& e) {
            std::cout << "Caught exception for invalid reshape: " << e.what() << "\n\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << std::endl;
    }

    return 0;
}

int test1() {
    try {
        // 1. Create a 4D array with shape (2, 3, 5, 4)
        std::vector<int> data(2 * 3 * 5 * 4);
        // Initialize data with sequential values for demonstration
        std::iota(data.begin(), data.end(), 0);

        MemoryView<int> mem1(data.data(), {2, 3, 5, 4});
        std::cout << "Original shape: (2, 3, 5, 4)" << std::endl;

        // 2. Perform slicing: mem1.slice({{}, {1}, {1, 3}, {}})
        // Equivalent to mem1[:, 1, 1:3, :]
        std::vector<Slice> slices = { Slice(), Slice(1), Slice(1, 3), Slice() };
        MemoryView<int> mem3 = mem1.slice(slices);
        std::cout << "After slicing with {{}, {1}, {1, 3}, {}}:" << std::endl;
        std::vector<size_t> shape3 = mem3.shape();
        std::cout << "Shape: (";
        for (size_t i = 0; i < shape3.size(); ++i) {
            std::cout << shape3[i];
            if (i != shape3.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // 3. Reshape mem3 to (2, 2, 4, 1)
        mem3.reshape({2, 2, 4, 1});
        std::cout << "After reshaping mem3 to (2, 2, 4, 1):" << std::endl;
        std::vector<size_t> shape_reshaped = mem3.shape();
        std::cout << "Shape: (";
        for (size_t i = 0; i < shape_reshaped.size(); ++i) {
            std::cout << shape_reshaped[i];
            if (i != shape_reshaped.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // 4. Access an element: mem3( {1, 1, 3, 0} )
        int value = mem3( {1, 1, 3, 0} );
        std::cout << "Element at (1, 1, 3, 0): " << value << std::endl;

        // 5. Slice with negative indices and partial slicing
        // mem1.slice({-1, {}, {}, {2}})
        std::vector<Slice> slices_neg = { Slice(-1), Slice(), Slice(), Slice(2) };
        MemoryView<int> mem_neg = mem1.slice(slices_neg);
        std::cout << "After slicing with {-1, {}, {}, {2}}:" << std::endl;
        std::vector<size_t> shape_neg = mem_neg.shape();
        std::cout << "Shape: (";
        for (size_t i = 0; i < shape_neg.size(); ++i) {
            std::cout << shape_neg[i];
            if (i != shape_neg.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // 6. Overload [] operator
        // mem1[1] corresponds to mem1.slice({1, {}, {}, {}})
        MemoryView<int> mem2 = mem1[1];
        std::cout << "After mem1[1]:" << std::endl;
        std::vector<size_t> shape2 = mem2.shape();
        std::cout << "Shape: (";
        for (size_t i = 0; i < shape2.size(); ++i) {
            std::cout << shape2[i];
            if (i != shape2.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // 7. Attempt an invalid slice (out of bounds)
        try {
            std::vector<Slice> invalid_slices = { Slice(0), Slice(5) }; // Second dimension size is 3
            MemoryView<int> mem_invalid = mem1.slice(invalid_slices);
        } catch (const std::exception& e) {
            std::cout << "Caught exception for invalid slice: " << e.what() << std::endl;
        }

        // 8. Attempt an invalid reshape
        try {
            mem1.reshape({4, 4}); // Total elements 2*3*5*4 = 120 vs 4*4 = 16
        } catch (const std::exception& e) {
            std::cout << "Caught exception for invalid reshape: " << e.what() << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << std::endl;
    }

    return 0;
}
int main()
{

    test1();
    test2();
}