#pragma once 

#include "../camera/camera.hpp"


/* ============== Plot DK_hierarchy animation ===============*/
struct PointKey {
    double x,y,z;
    bool operator<(const PointKey& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
};

struct VertexInfo {
    int id;
    Point p;
};

struct FaceInfo {
    int id;
    std::array<int,3> vertices;
};

PointKey key(const Point& p);

void extract_polyhedron(
    const Polyhedron& poly,
    std::map<PointKey,int>& point_id,
    std::vector<VertexInfo>& vertices,
    std::vector<FaceInfo>& faces,
    int& next_vid,
    int& next_fid);


/*========= Plot intersection of polyhedra ================*/
void polyhedron_to_json(
    const Polyhedron& poly,
    json& Jpoly
);

void export_scene_json(
    const Polyhedron& P,
    const Polyhedron& Q,
    const Point& O,
    const bool intersects,
    const std::string& filename
);

/* same as the above function, adding color intersection*/
void export_intersection_json(
    std::vector<Polyhedron>& fov_polys,
    const Polyhedron& target_intersection,
    const Point& O,
    const std::string& filename
) ;
/* ====== Plot Chazelle's definition for supprting alg =======*/

void export_dual_points(
    const std::vector<Point>& pts,
    const std::string& file,
    const std::string& color
) ;


