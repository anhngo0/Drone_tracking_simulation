#include "../../include/camera/convex_polyhedron.hpp"

// triangulation
void triangulate_in_place(Polyhedron& poly) {
    CGAL::Polygon_mesh_processing::triangulate_faces(poly);
}

bool is_point_inside_polyhedron(
    const Polyhedron& poly,
    const Point& query)
{
    if (poly.is_empty())
        return false;

    // Lấy một điểm bên trong để xác định orientation
    Point centroid = CGAL::centroid(
        poly.points_begin(),
        poly.points_end());

    for (auto f = poly.facets_begin();
         f != poly.facets_end(); ++f)
    {
        auto h = f->halfedge();

        Point p = h->vertex()->point();
        Point q = h->next()->vertex()->point();
        Point r = h->next()->next()->vertex()->point();

        Plane pl(p, q, r);

        // Đảm bảo negative side là inside
        if (pl.oriented_side(centroid)
            != CGAL::ON_NEGATIVE_SIDE)
        {
            pl = pl.opposite();
        }

        // Nếu điểm nằm phía positive side → ngoài
        if (pl.oriented_side(query)
            == CGAL::ON_POSITIVE_SIDE)
        {
            return false;
        }
    }

    return true; // nằm trong hoặc trên biên
}
/*==============================================================*/
/*=== CHECKING CONTAINMENT USING BOUNDING BOX (FAST CHECK) =====*/
/*==============================================================*/

AABB compute_aabb(const Polyhedron& P) {
    AABB box;

    auto it = P.points_begin();
    box.xmin = box.xmax = it->x();
    box.ymin = box.ymax = it->y();
    box.zmin = box.zmax = it->z();

    ++it;
    for (; it != P.points_end(); ++it) {
        box.xmin = std::min(box.xmin, it->x());
        box.xmax = std::max(box.xmax, it->x());
        box.ymin = std::min(box.ymin, it->y());
        box.ymax = std::max(box.ymax, it->y());
        box.zmin = std::min(box.zmin, it->z());
        box.zmax = std::max(box.zmax, it->z());
    }
    return box;
}

bool box_inside_box(const AABB& P, const AABB& Q) {
    return
        P.xmin >= Q.xmin && P.xmax <= Q.xmax &&
        P.ymin >= Q.ymin && P.ymax <= Q.ymax &&
        P.zmin >= Q.zmin && P.zmax <= Q.zmax;
}

/*========================================================*/
/*=============  CHEKC CONTAINMENT  ============*/
/*========================================================*/

bool point_inside_face(const Point& p, const Plane& pl) {
    return pl.oriented_side(p) != CGAL::ON_POSITIVE_SIDE;
}

std::optional<Point> find_chebyshev_center_LP(
    std::vector<Polyhedron>& polys)
{

    // variables: x(0), y(1), z(2), r(3)
    Program lp(CGAL::SMALLER, false, 0, false, 0);

    // r >= 0
    lp.set_l(3, true, ET(0));

    int constraint_idx = 0;

    for (auto& poly : polys)
    {
        if (poly.is_empty())
            return std::nullopt;

        // centroid để xác định orientation
        Point centroid = find_volume_centroid(poly);

        for (auto f = poly.facets_begin();
             f != poly.facets_end(); ++f)
        {
            auto h = f->halfedge();

            Point p = h->vertex()->point();
            Point q = h->next()->vertex()->point();
            Point r_pt = h->next()->next()->vertex()->point();

            Plane pl(p, q, r_pt);

            // đảm bảo negative side là inside
            if (pl.oriented_side(centroid)
                != CGAL::ON_NEGATIVE_SIDE)
            {
                pl = pl.opposite();
            }

            double a = CGAL::to_double(pl.a());
            double b = CGAL::to_double(pl.b());
            double c = CGAL::to_double(pl.c());
            double d = CGAL::to_double(pl.d());

            // ||n||
            ET norm = CGAL::sqrt(a*a + b*b + c*c);

            // a x + b y + c z + r*norm <= -d
            lp.set_a(0, constraint_idx, a);
            lp.set_a(1, constraint_idx, b);
            lp.set_a(2, constraint_idx, c);
            lp.set_a(3, constraint_idx, norm);

            lp.set_b(constraint_idx, -d);

            constraint_idx++;
        }
    }

    // maximize r  <=> minimize -r
    lp.set_c(3, ET(-1));

    Solution s = CGAL::solve_linear_program(lp, ET());

    if (!s.is_optimal())
        return std::nullopt;

    auto it = s.variable_values_begin();

    double x = CGAL::to_double(*it++);
    double y = CGAL::to_double(*it++);
    double z = CGAL::to_double(*it++);
    double r = CGAL::to_double(*it);

    if (r < ET(0))   // cho phép r = 0 nếu muốn giữ giao suy biến
        return std::nullopt;

    return Point(x, y, z);
}

/*========================================================*/
/*========  COMPUTE INTERSECTION OF POLYHEDRA  ===========*/
/*========================================================*/

Polyhedron compute_intersection_LP(std::vector<Polyhedron>& polys)
{

    // 1️⃣ Tìm Chebyshev center (interior point)
    auto center_opt = find_chebyshev_center_LP(polys);

    if (!center_opt.has_value())
        return Polyhedron();   // không giao hoặc giao suy biến

    Point interior_point = center_opt.value();

    // 2️⃣ Thu thập toàn bộ halfspaces từ P và Q
    std::vector<Plane> planes;

    auto extract_planes = [&](Polyhedron& poly)
    {
        // centroid để xác định orientation
        Point centroid = find_volume_centroid(poly);

        for (auto f = poly.facets_begin();
             f != poly.facets_end(); ++f)
        {
            auto h = f->halfedge();

            Point p = h->vertex()->point();
            Point q = h->next()->vertex()->point();
            Point r = h->next()->next()->vertex()->point();

            Plane pl(p, q, r);

            // đảm bảo negative side là inside
            if (pl.oriented_side(centroid)
                != CGAL::ON_NEGATIVE_SIDE)
            {
                pl = pl.opposite();
            }

            planes.push_back(pl);
        }
    };

    for(Polyhedron poly_fov: polys){
        extract_planes(poly_fov);
    }

    // 3️⃣ Đảm bảo interior_point nằm ở negative side
    for (auto& pl : planes)
    {
        if (pl.oriented_side(interior_point)
            != CGAL::ON_NEGATIVE_SIDE &&
            pl.oriented_side(interior_point)
            != CGAL::ON_BOUNDARY)
        {
            pl = pl.opposite();
        }
    }

    // 4️⃣ Dựng polyhedron giao bằng halfspace intersection
    Polyhedron result;

    CGAL::halfspace_intersection_3(
        planes.begin(),
        planes.end(),
        result,
        interior_point);

    if (!result.is_empty())
        CGAL::Polygon_mesh_processing::triangulate_faces(result);

    return result;
}

/*=================================================================================*/
/*============= INTERSECTION OF POLYHEDRA ============*/
/*=================================================================================*/

std::vector<Tetrahedron> decompose_to_tetrahedra(
    Polyhedron& poly, 
    const Point& M) 
{
    triangulate_in_place(poly);

    std::vector<Tetrahedron> tets;

    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        // get 3 vertices of a triangle
        auto h = f->halfedge();
        Point p1 = h->vertex()->point();
        Point p2 = h->next()->vertex()->point();
        Point p3 = h->next()->next()->vertex()->point();

        // create a polyhedron from 3 vertices and an inner point
        double vol = CGAL::to_double(CGAL::volume(M, p1, p2, p3));
        
        if (std::abs(vol) > 1e-14) {
            tets.push_back({M, p1, p2, p3, std::abs(vol)});
        }
    }
    return tets;
}

Point find_volume_centroid(Polyhedron& poly) {

    double total_vol = 0.0;
    double cx = 0, cy = 0, cz = 0;
    
    // Choose a fixed reference point (e.g., the first vertex)
    Point origin = poly.vertices_begin()->point();

    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        auto h = f->halfedge();
        Point p1 = h->vertex()->point();
        Point p2 = h->next()->vertex()->point();
        Point p3 = h->next()->next()->vertex()->point();

        // Signed volume of the tetrahedron formed by origin and the triangular face
        double v = CGAL::to_double(CGAL::volume(origin, p1, p2, p3));
        
        // Centroid of the tetrahedron: (p1+p2+p3+origin)/4
        total_vol += v;
        cx += v * (CGAL::to_double(p1.x() + p2.x() + p3.x() + origin.x()) / 4.0);
        cy += v * (CGAL::to_double(p1.y() + p2.y() + p3.y() + origin.y()) / 4.0);
        cz += v * (CGAL::to_double(p1.z() + p2.z() + p3.z() + origin.z()) / 4.0);
    }

    if (std::abs(total_vol) < 1e-12) return origin; 

    // Volume-weighted centroid of the polyhedron
    return Point(cx / total_vol, cy / total_vol, cz / total_vol);
}

// Helper: convert degrees to radians
inline double deg2rad(double deg) { return deg * M_PI / 180.0; }


