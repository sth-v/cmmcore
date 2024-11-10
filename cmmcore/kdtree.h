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

#include "vec.h"
namespace cmmcore{

// Structure to represent a 3D point


// KD-Tree Node structure
template<typename PT>
struct KDNode {
    PT point;
    //size_t index; // Index of the point in the original point list
    int axis; // 0: x, 1: y, 2: z
    KDNode* left;
    KDNode* right;

    KDNode(const PT& pt, int ax) : point(pt), axis(ax), left(nullptr), right(nullptr) {}
};

// Comparator for sorting points based on a specific axis
struct PointComparator3D {
    int axis;
    PointComparator3D(int ax) : axis(ax) {}
    bool operator()(const vec3& a, const vec3& b) const {
        if (axis == 0) return a.x < b.x;
        if (axis == 1) return a.y < b.y;
        return a.z < b.z;
    }
};

// KD-Tree Class
template<typename PT, typename Comparator>
class KDTree {
    using node_type=KDNode<PT>;
private:
    node_type* root;
    size_t size_;

    // Recursive function to build KD-Tree
    node_type* buildKDTree(typename std::vector<PT>::iterator begin, typename std::vector<PT>::iterator end, int depth) {
        if (begin >= end)
            return nullptr;

        int axis = depth % 3;
        size_t len = distance(begin, end);
        auto mid = begin + len / 2;

        // Partially sort the points to find the median on current axis
        std::nth_element(begin, mid, end, Comparator(axis));

        // Create node and construct subtrees
        node_type* node = new KDNode(*mid, axis);
        node->left = buildKDTree(begin, mid, depth + 1);
        node->right = buildKDTree(mid + 1, end, depth + 1);
        return node;
    }

    // Recursive function to insert a point into the KD-Tree
    node_type* insertRec(node_type* node, const PT& point, int depth) {
        if (node == nullptr) {
            size_++;
            return new node_type(point, depth % 3);
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
    void rangeSearchRec(node_type* node, const PT& lower, const PT& upper, std::vector<PT>& result) const {
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
    void nearestNeighborRec(node_type* node, const PT& target, PT& best, double& bestDistSq) const {
        if (node == nullptr)
            return;

        double d = target.distanceSq(node->point);
        if (d < bestDistSq) {
            bestDistSq = d;
            best = node->point;
        }

        int axis = node->axis;
        node_type* nearChild = nullptr;
        node_type* farChild = nullptr;

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
    void deleteKDTree(node_type* node) {
        if (node == nullptr)
            return;
        deleteKDTree(node->left);
        deleteKDTree(node->right);
        delete node;
    }

public:
    // Constructor
    KDTree() : root(nullptr), size_(0) {}

    // Destructor
    ~KDTree() {
        deleteKDTree(root);
    }

    // Build the KD-Tree from a list of points
    void build(std::vector<PT>& points) {
        deleteKDTree(root); // Clean existing tree
        size_ = points.size();
        root = buildKDTree(points.begin(), points.end(), 0);
    }

    // Insert a new point into the KD-Tree
    void insert(const PT& point) {
        root = insertRec(root, point, 0);
    }

    // Perform a range search within the bounding box defined by lower and upper points
    std::vector<PT> rangeSearch(const PT& lower, const PT& upper) const {
        std::vector<PT> result;
        rangeSearchRec(root, lower, upper, result);
        return result;
    }

    // Find the nearest neighbor to the target point
    PT nearestNeighbor(const PT& target) const {
        if (root == nullptr)
            throw std::runtime_error("KD-Tree is empty!");

        PT best = root->point;
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

// Utility function to print a PT
void printPoint(const vec3& p) {
    std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}
    typedef KDTree<vec3,PointComparator3D> KDTree3D;

}
#endif //CMMCORE_KDTREE_H
