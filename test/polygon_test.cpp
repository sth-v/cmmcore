//
// Created by Andrew Astakhov on 08.10.24.
//
#include <iostream>
#include <vector>
#include "cmmcore/polygon.h"
using namespace cmmcore;
// Assume the previously defined code (vec2, Vector, PolygonRelationship2D enum, and polygonPolygonRelationship2D function) is included here

void printPolygonRelationship2D(PolygonRelationship2D rel) {
    switch (rel) {
        case PolygonRelationship2D::INTERSECT:
            std::cout << "The polygons intersect." << std::endl;
        break;
        case PolygonRelationship2D::TOUCH:
            std::cout << "The polygons touch." << std::endl;
        break;
        case PolygonRelationship2D::DISTANCE:
            std::cout << "The polygons are separate." << std::endl;
        break;
    }
}

int main() {
    // Example 1: Intersecting polygons
    std::vector<vec2> polygon1 = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    std::vector<vec2> polygon2 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
    
    std::cout << "Example 1:" << std::endl;
    printPolygonRelationship2D(polygonRelationship2D(polygon1, polygon2));

    // Example 2: Touching polygons
    std::vector<vec2> polygon3 = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    std::vector<vec2> polygon4 = {{2, 0}, {4, 0}, {4, 2}, {2, 2}};
    
    std::cout << "\nExample 2:" << std::endl;
    printPolygonRelationship2D(polygonRelationship2D(polygon3, polygon4));

    // Example 3: Separate polygons
    std::vector<vec2> polygon5 = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    std::vector<vec2> polygon6 = {{3, 3}, {5, 3}, {5, 5}, {3, 5}};
    
    std::cout << "\nExample 3:" << std::endl;
    printPolygonRelationship2D(polygonRelationship2D(polygon5, polygon6));

    // Example 4: One polygon inside another
    std::vector<vec2> polygon7 = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    std::vector<vec2> polygon8 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
    
    std::cout << "\nExample 4:" << std::endl;
    printPolygonRelationship2D(polygonRelationship2D(polygon7, polygon8));

    return 0;
}