//
// Created by Andrew Astakhov on 05.11.24.
//

#ifndef CMMCORE_KDTREE_H
#define CMMCORE_KDTREE_H
#include <iostream>
#include <array>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>
#include "vec.h"
namespace cmmcore{

// Structure to represent a 3D point


// KD-Tree Node structure
struct KDNode {
    vec3 point;
    int axis; // 0: x, 1: y, 2: z
    KDNode* left;
    KDNode* right;

    KDNode(const vec3& pt, int ax) : point(pt), axis(ax), left(nullptr), right(nullptr) {}
};

// Comparator for sorting points based on a specific axis
struct PointComparator {
    int axis;
    PointComparator(int ax) : axis(ax) {}
    bool operator()(const vec3& a, const vec3& b) const {
        if (axis == 0) return a.x < b.x;
        if (axis == 1) return a.y < b.y;
        return a.z < b.z;
    }
};

// KD-Tree Class
class KDTree3D {
private:
    KDNode* root;
    size_t size_;

    // Recursive function to build KD-Tree
    KDNode* buildKDTree(std::vector<vec3>::iterator begin, std::vector<vec3>::iterator end, int depth) {
        if (begin >= end)
            return nullptr;

        int axis = depth % 3;
        size_t len = distance(begin, end);
        auto mid = begin + len / 2;

        // Partially sort the points to find the median on current axis
        nth_element(begin, mid, end, PointComparator(axis));

        // Create node and construct subtrees
        KDNode* node = new KDNode(*mid, axis);
        node->left = buildKDTree(begin, mid, depth + 1);
        node->right = buildKDTree(mid + 1, end, depth + 1);
        return node;
    }

    // Recursive function to insert a point into the KD-Tree
    KDNode* insertRec(KDNode* node, const vec3& point, int depth) {
        if (node == nullptr) {
            size_++;
            return new KDNode(point, depth % 3);
        }

        int axis = node->axis;
        if ((axis == 0 && point.x < node->point.x) ||
            (axis == 1 && point.y < node->point.y) ||
            (axis == 2 && point.z < node->point.z)) {
            node->left = insertRec(node->left, point, depth + 1);
        }
        else {
            node->right = insertRec(node->right, point, depth + 1);
        }
        return node;
    }

    // Recursive function for range search
    void rangeSearchRec(KDNode* node, const vec3& lower, const vec3& upper, std::vector<vec3>& result) const {
        if (node == nullptr)
            return;

        // Check if current node is inside the range
        if (node->point.x >= lower.x && node->point.x <= upper.x &&
            node->point.y >= lower.y && node->point.y <= upper.y &&
            node->point.z >= lower.z && node->point.z <= upper.z) {
            result.push_back(node->point);
        }

        // Decide whether to traverse left/right subtree
        int axis = node->axis;
        if ((axis == 0 && lower.x <= node->point.x) ||
            (axis == 1 && lower.y <= node->point.y) ||
            (axis == 2 && lower.z <= node->point.z)) {
            rangeSearchRec(node->left, lower, upper, result);
        }
        if ((axis == 0 && upper.x >= node->point.x) ||
            (axis == 1 && upper.y >= node->point.y) ||
            (axis == 2 && upper.z >= node->point.z)) {
            rangeSearchRec(node->right, lower, upper, result);
        }
    }

    // Recursive function for nearest neighbor search
    void nearestNeighborRec(KDNode* node, const vec3& target, vec3& best, double& bestDistSq) const {
        if (node == nullptr)
            return;

        double d = target.distanceSq(node->point);
        if (d < bestDistSq) {
            bestDistSq = d;
            best = node->point;
        }

        int axis = node->axis;
        KDNode* nearChild = nullptr;
        KDNode* farChild = nullptr;

        if ((axis == 0 && target.x < node->point.x) ||
            (axis == 1 && target.y < node->point.y) ||
            (axis == 2 && target.z < node->point.z)) {
            nearChild = node->left;
            farChild = node->right;
        }
        else {
            nearChild = node->right;
            farChild = node->left;
        }

        nearestNeighborRec(nearChild, target, best, bestDistSq);

        // Check if we need to explore the far child
        double diff = 0.0;
        if (axis == 0) diff = target.x - node->point.x;
        else if (axis == 1) diff = target.y - node->point.y;
        else diff = target.z - node->point.z;

        if (diff * diff < bestDistSq) {
            nearestNeighborRec(farChild, target, best, bestDistSq);
        }
    }

    // Recursive function to delete the KD-Tree
    void deleteKDTree(KDNode* node) {
        if (node == nullptr)
            return;
        deleteKDTree(node->left);
        deleteKDTree(node->right);
        delete node;
    }

public:
    // Constructor
    KDTree3D() : root(nullptr), size_(0) {}

    // Destructor
    ~KDTree3D() {
        deleteKDTree(root);
    }

    // Build the KD-Tree from a list of points
    void build(std::vector<vec3>& points) {
        deleteKDTree(root); // Clean existing tree
        size_ = points.size();
        root = buildKDTree(points.begin(), points.end(), 0);
    }

    // Insert a new point into the KD-Tree
    void insert(const vec3& point) {
        root = insertRec(root, point, 0);
    }

    // Perform a range search within the bounding box defined by lower and upper points
    std::vector<vec3> rangeSearch(const vec3& lower, const vec3& upper) const {
        std::vector<vec3> result;
        rangeSearchRec(root, lower, upper, result);
        return result;
    }

    // Find the nearest neighbor to the target point
    vec3 nearestNeighbor(const vec3& target) const {
        if (root == nullptr)
            throw std::runtime_error("KD-Tree is empty!");

        vec3 best = root->point;
        double bestDistSq = target.distanceSq(best);
        nearestNeighborRec(root, target, best, bestDistSq);
        return best;
    }

    // Get the number of points in the KD-Tree
    size_t size() const {
        return size_;
    }

    // Check if the KD-Tree is empty
    bool empty() const {
        return size_ == 0;
    }
};

// Utility function to print a vec3
void printPoint(const vec3& p) {
    std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}


}
#endif //CMMCORE_KDTREE_H
