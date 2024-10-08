/*
Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef POLYLINE_H
#define POLYLINE_H
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <cmath>
#include <limits>
#ifdef CYTHON_ABI
#include "vec.h"
#else
#include "cmmcore/vec.h"
#endif

namespace cmmcore{

// Define epsilon for distance comparison
const double eps = std::numeric_limits<double>::epsilon();

// Quantization function for spatial hashing
struct PointKey {
    int x, y, z;
    PointKey() = default;
    PointKey(double x, double y, double z) : x(x), y(y), z(z) {}


      PointKey(const vec3& point) {
        x = static_cast<int>(std::floor(point.x / eps + 0.5));
        y = static_cast<int>(std::floor(point.y / eps + 0.5));
        z = static_cast<int>(std::floor(point.z / eps + 0.5));
    }


    bool operator==(const PointKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Custom hash function for PointKey
struct PointKeyHash {
    std::size_t operator()(const PointKey& key) const {
        return (static_cast<std::size_t>(key.x) * 73856093) ^
               (static_cast<std::size_t>(key.y) * 19349663) ^
               (static_cast<std::size_t>(key.z) * 83492791);
    }
};

// Segment structure
struct Segment {
    vec3 p1;
    vec3 p2;
};

// Endpoint structure
struct Endpoint {
    vec3 point;
    Segment* segment;
    bool is_p1; // true if this endpoint is p1 of the segment, false if p2
};

// Spatial index mapping quantized points to endpoints
using SpatialIndex = std::unordered_map<PointKey, std::vector<Endpoint*>, PointKeyHash>;

// Function to add an endpoint to the spatial index
void add_to_index(SpatialIndex& index, Endpoint* endpoint) {
    PointKey key(endpoint->point);
    index[key].push_back(endpoint);
}

// Function to remove an endpoint from the spatial index
void remove_from_index(SpatialIndex& index, Endpoint* endpoint) {
    PointKey key={endpoint->point.x, endpoint->point.y, endpoint->point.z};
    auto& vec = index[key];
    vec.erase(std::remove(vec.begin(), vec.end(), endpoint), vec.end());
    if (vec.empty()) {
        index.erase(key);
    }
}

// Function to find endpoints near a given point
std::vector<Endpoint*> find_nearby_endpoints(const SpatialIndex& index, const vec3& point) {
    std::vector<Endpoint*> result;
    PointKey key(point);
    // Check neighboring keys (27 cells in total)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                PointKey neighbor_key(key.x + dx, key.y + dy, key.z + dz);
                auto it = index.find(neighbor_key);
                if (it != index.end()) {
                    result.insert(result.end(), it->second.begin(), it->second.end());
                }
            }
        }
    }
    return result;
}

// Main function to assemble polylines from segments
 std::vector<std::list<vec3>>& assemble_polylines(std::vector<Segment>& segments,std::vector<std::list<vec3>>& branches) {
    // Create spatial index and endpoint list
    SpatialIndex index;
    std::vector<Endpoint> endpoints;

    // Initialize endpoints and spatial index
    for (auto& segment : segments) {
        endpoints.push_back({segment.p1, &segment, true});
        endpoints.push_back({segment.p2, &segment, false});
    }

    for (auto& endpoint : endpoints) {
        add_to_index(index, &endpoint);
    }

    // Polylines result


    // Set to keep track of processed segments
    std::unordered_set<Segment*> processed_segments;

    // Main loop
    for (auto& endpoint : endpoints) {
        if (processed_segments.count(endpoint.segment)) {
            continue; // Skip if segment is already processed
        }

        // Start a new polyline
        std::list<vec3> polyline;
        polyline.push_back(endpoint.point);
        polyline.push_back(endpoint.is_p1 ? endpoint.segment->p2 : endpoint.segment->p1);

        processed_segments.insert(endpoint.segment);
        remove_from_index(index, &endpoint);
        // Remove the opposite endpoint
        Endpoint* opposite_endpoint = nullptr;
        for (auto& ep : endpoints) {
            if (ep.segment == endpoint.segment && ep.is_p1 != endpoint.is_p1) {
                opposite_endpoint = &ep;
                break;
            }
        }
        if (opposite_endpoint) {
            remove_from_index(index, opposite_endpoint);
        }

        // Extend polyline from both ends
        bool extended = true;
        while (extended) {
            extended = false;
            // Try to extend at the front
            vec3& current_point_front = polyline.front();
            auto nearby_endpoints_front = find_nearby_endpoints(index, current_point_front);
            for (auto* ep : nearby_endpoints_front) {
                if (processed_segments.count(ep->segment)) {
                    continue;
                }
                if ((ep->point - current_point_front).length() < eps) {
                    // Found a connected segment
                    vec3 far_point = ep->is_p1 ? ep->segment->p2 : ep->segment->p1;
                    polyline.push_front(far_point);
                    processed_segments.insert(ep->segment);

                    // Remove endpoints from index
                    remove_from_index(index, ep);
                    // Remove the opposite endpoint
                    Endpoint* opp_ep = nullptr;
                    for (auto& e : endpoints) {
                        if (e.segment == ep->segment && e.is_p1 != ep->is_p1) {
                            opp_ep = &e;
                            break;
                        }
                    }
                    if (opp_ep) {
                        remove_from_index(index, opp_ep);
                    }

                    extended = true;
                    break; // Break to restart from the new front
                }
            }

            // Try to extend at the back
            vec3& current_point_back = polyline.back();
            auto nearby_endpoints_back = find_nearby_endpoints(index, current_point_back);
            for (auto* ep : nearby_endpoints_back) {
                if (processed_segments.count(ep->segment)) {
                    continue;
                }
                if ((ep->point - current_point_back).length() < eps) {
                    // Found a connected segment
                    vec3 far_point = ep->is_p1 ? ep->segment->p2 : ep->segment->p1;
                    polyline.push_back(far_point);
                    processed_segments.insert(ep->segment);

                    // Remove endpoints from index
                    remove_from_index(index, ep);
                    // Remove the opposite endpoint
                    Endpoint* opp_ep = nullptr;
                    for (auto& e : endpoints) {
                        if (e.segment == ep->segment && e.is_p1 != ep->is_p1) {
                            opp_ep = &e;
                            break;
                        }
                    }
                    if (opp_ep) {
                        remove_from_index(index, opp_ep);
                    }

                    extended = true;
                    break; // Break to restart from the new back
                }
            }
        }

        // Check if polyline is closed
        if ((polyline.front() - polyline.back()).length() < eps) {
            // Remove the duplicate point at the end
            polyline.pop_back();
            // Optionally, mark as closed if needed
        }

        // Add the polyline to the result
        branches.push_back(std::move(polyline));
    }

    return branches;
}
/*
// Example usage
int main() {
    // Sample segments forming two polylines (one closed, one open)
    std::vector<Segment> segments = {
        {{0, 0, 0}, {1, 0, 0}},
        {{1, 0, 0}, {1, 1, 0}},
        {{1, 1, 0}, {0, 1, 0}},
        {{0, 1, 0}, {0, 0, 0}}, // Closed loop

        {{2, 0, 0}, {3, 0, 0}},
        {{3, 0, 0}, {3, 1, 0}},
        {{3, 1, 0}, {2, 1, 0}}, // Open polyline
    };

    auto polylines = assemble_polylines(segments);

    // Output the polylines
    int polyline_index = 1;
    for (const auto& polyline : polylines) {
        std::cout << "Polyline " << polyline_index++ << " (" << (polyline.front() == polyline.back() ? "Closed" : "Open") << "):" << std::endl;
        for (const auto& point : polyline) {
            std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ") ";
        }
        std::cout << std::endl;
    }

    return 0;
}
*/
}
#endif //POLYLINE_H
