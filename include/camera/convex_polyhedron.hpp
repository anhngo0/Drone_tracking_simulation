#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Gmpq.h>

#include <nlohmann/json.hpp>
#include <vector>
#include <set>
#include <random>
#include <algorithm>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point  = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Plane  = Kernel::Plane_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

using FT = Kernel::FT;
using ET = CGAL::Gmpq;
using Program = CGAL::Quadratic_program<ET>;
using Solution = CGAL::Quadratic_program_solution<ET>;

using json = nlohmann::json;

using ScalarFunction = std::function<double(double, double, double)>;

bool is_point_inside_polyhedron(
    const Polyhedron& poly,
    const Point& query);

Polyhedron triangulate_in_place(const Polyhedron& poly);

/*==============================================================*/
/*=== CHECKING CONTAINMENT USING BOUNDING BOX (FAST CHECK) =====*/
/*==============================================================*/

struct AABB {
    FT xmin, xmax;
    FT ymin, ymax;
    FT zmin, zmax;
};

/*==============================================================*/
/*============= CHECKING INTERSECTION OF POLYHEDRA =============*/
/*==============================================================*/

struct IntersectionResult {
    bool intersects;
    Point point;
};

/*========================================================*/
/*========  COMPUTE INTERSECTION OF POLYHEDRA  ===========*/
/*========================================================*/

Polyhedron compute_intersection_LP(std::vector<Polyhedron>& polys);
/*=================================================================================*/
/*=============  COMPUTE GAUSS ITERATION IN INTERSECTION OF POLYHEDRA ============*/
/*=================================================================================*/
struct Tetrahedron {
    Point v0, v1, v2, v3;
    double volume;
};

ScalarFunction make_gaussian_function(Point mean, double delta);

Point find_volume_centroid(Polyhedron& poly);

