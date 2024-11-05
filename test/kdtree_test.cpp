//
// Created by Andrew Astakhov on 05.11.24.
//
// Example usage

#include "cmmcore/kdtree.h"
using namespace cmmcore;
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Sample points
    std::vector<vec3> points = {
        {2.0, 3.0, 4.0},
        {5.0, 4.0, 2.0},
        {9.0, 6.0, 7.0},
        {4.0, 7.0, 9.0},
        {8.0, 1.0, 5.0},
        {7.0, 2.0, 6.0}
    };

    // Initialize KD-Tree and build it
    KDTree3D tree;
    tree.build(points);
    std::cout << "KD-Tree built with " << tree.size() << " points.\n";

    // Insert a new point
    vec3 newPoint(3.0, 5.0, 4.0);
    tree.insert(newPoint);
    std::cout << "Inserted point ";
    printPoint(newPoint);
    std::cout << ". KD-Tree now has " << tree.size() << " points.\n";

    // Perform a range search
    vec3 lower(2.0, 2.0, 2.0);
    vec3 upper(6.0, 6.0, 6.0);
    std::vector<vec3> rangeResult = tree.rangeSearch(lower, upper);
    std::cout << "Range Search between ";
    printPoint(lower);
    std::cout << " and ";
    printPoint(upper);
    std::cout << " found " << rangeResult.size() << " points:\n";
    for (const auto& p : rangeResult) {
        printPoint(p);
        std::cout << "\n";
    }

    // Perform nearest neighbor search
    vec3 target(9.0, 2.0, 6.0);
    try {
        vec3 nearest = tree.nearestNeighbor(target);
        std::cout << "Nearest neighbor to ";
        printPoint(target);
        std::cout << " is ";
        printPoint(nearest);
        std::cout << " with squared distance " << target.distanceSq(nearest) << ".\n";
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
    }

    return 0;
}