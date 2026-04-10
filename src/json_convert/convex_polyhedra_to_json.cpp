#include "../../include/json_convert/convex_polyhedra_to_json.hpp"


/*plot DK hierarchy animation*/
PointKey key(const Point& p) {
    double point_x = CGAL::to_double(p.x());
    double point_y = CGAL::to_double(p.y());
    double point_z = CGAL::to_double(p.z());
    return {point_x,point_y,point_z};
}

void extract_polyhedron(
    const Polyhedron& poly,
    std::map<PointKey,int>& point_id,
    std::vector<VertexInfo>& vertices,
    std::vector<FaceInfo>& faces,
    int& next_vid,
    int& next_fid)
{
    // ===== vertices =====
    for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
        PointKey k = key(v->point());
        if (!point_id.count(k)) {
            point_id[k] = next_vid;
            vertices.push_back({
                next_vid,
                v->point()
            });
            ++next_vid;
        }
    }

    // ===== faces (poly MUST be triangulated) =====
    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {

        auto h = f->facet_begin();
        std::array<int,3> vids;
        int i = 0;

        do {
            PointKey k = key(h->vertex()->point());
            vids[i++] = point_id[k];
        } while (++h != f->facet_begin() && i < 3);

        // chỉ nhận tam giác
        if (i == 3) {
            faces.push_back({
                next_fid++,
                vids
            });
        }
    }
}

/*Plot intersection of polyhedra*/

void polyhedron_to_json(
    const Polyhedron& poly,
    json& Jpoly
) {
    std::map<Polyhedron::Vertex_const_handle, int> vid;
    int vcount = 0, fcount = 0;

    Jpoly["vertices"] = json::array();
    Jpoly["faces"]    = json::array();

    // ===== vertices =====
    for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
        const Point& p = v->point();
        vid[v] = vcount;

        Jpoly["vertices"].push_back({
            {"id", vcount},
            {"pos", {
                CGAL::to_double(p.x()),
                CGAL::to_double(p.y()),
                CGAL::to_double(p.z())
            }}
        });

        ++vcount;
    }

    // ===== faces (assumed triangulated) =====
    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        auto h = f->facet_begin();

        std::array<int,3> fv;
        int i = 0;
        do {
            fv[i++] = vid[h->vertex()];
        } while (++h != f->facet_begin() && i < 3);

        if (i == 3) {
            Jpoly["faces"].push_back({
                {"id", fcount++},
                {"v", { fv[0], fv[1], fv[2] }}
            });
        }
    }
}

/*Plot 2 polyhedra with a point inside their intersection*/
void export_scene_json(
    const Polyhedron& P,
    const Polyhedron& Q,
    const Point& O,
    const bool intersects,
    const std::string& filename
) {
    json J;

    // ===== polyhedra =====
    J["polyhedra"] = json::object();

    polyhedron_to_json(P, J["polyhedra"]["P"]);
    polyhedron_to_json(Q, J["polyhedra"]["Q"]);

    // ===== intersection =====

    J["intersection"]["exists"] = intersects;

    if (intersects) {
        J["intersection"]["O"] = { 
            CGAL::to_double(O.x()), 
            CGAL::to_double(O.y()), 
            CGAL::to_double(O.z()) 
        };
    }

    std::ofstream(filename) << J.dump(2);
}

