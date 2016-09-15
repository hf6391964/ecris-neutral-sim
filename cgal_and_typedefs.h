#ifndef CGAL_AND_TYPEDEFS_H
#define CGAL_AND_TYPEDEFS_H

#include <random>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include "STL_reader.h"
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Direction_3 Direction;
typedef K::Ray_3 Ray;
typedef CGAL::Bbox_3 Bbox;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef Surface_mesh::Face_index face_descriptor;
typedef Surface_mesh::Vertex_index vertex_descriptor;
typedef Surface_mesh::Halfedge_index halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type>
    Ray_intersection;

typedef std::mt19937 Rng;
static std::uniform_real_distribution<double> uni01(0.0, 1.0);

#endif

