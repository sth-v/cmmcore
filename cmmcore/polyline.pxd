cdef extern from "polyline.h"  namespace "cmmcore" nogil:
    
    # Forward declarations
    cdef cppclass vec3:
        pass
    
    cdef cppclass SpatialIndex:
        pass

    cdef cppclass list[T]:
        pass
        
    cdef cppclass vector[T]:
        void push_back(T&)

    # PointKey Structure
    cdef cppclass PointKey:
        int x, y, z
        PointKey() except +
        PointKey(double x, double y, double z) except +
        PointKey(const vec3& point) except +
        bool operator==(const PointKey& other) const

    # PointKeyHash Structure
    cdef cppclass PointKeyHash:
        size_t operator()(const PointKey& key) const

    # Segment Structure
    cdef cppclass Segment:
        vec3 p1
        vec3 p2

    # Endpoint Structure
    cdef cppclass Endpoint:
        vec3 point
        Segment* segment
        bool is_p1

    # Functions
    void add_to_index(SpatialIndex& index, Endpoint* endpoint)
    void remove_from_index(SpatialIndex& index, Endpoint* endpoint)
    vector[Endpoint*] find_nearby_endpoints(const SpatialIndex& index, const vec3& point)
    vector[list[vec3]]& assemble_polylines(vector[Segment]& segments, vector[list[vec3]]& branches)