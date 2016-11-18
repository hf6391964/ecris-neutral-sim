#pragma once

#include <random>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Direction_3 Direction;
typedef K::Ray_3 Ray;
typedef K::Aff_transformation_3 Aff_transformation;
typedef CGAL::Bbox_3 Bbox;

typedef std::mt19937 Rng;
static std::uniform_real_distribution<double> uni01(0.0, 1.0);

