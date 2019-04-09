
#pragma once

#include <Eigen/Dense>

#include "ArrangementDefs.h"
#include "PlanarArrangement.h"

class PatternBoundary {
public:
	PatternBoundary(const std::vector<Polygon_2>& boundary_polygons);
	PatternBoundary(const std::vector<Polyline_2>& boundary_polylines);
	Polyline_2 filter_and_snap(Polyline_2& polyline);
private:
	Polyline_2 filter(Polyline_2& polyline);
	Polyline_2 snap(Polyline_2& polyline);
	bool inside_boundary(const Point_2& pt);

	static std::vector<Polygon_2> polylines_to_polygons(const std::vector<Polyline_2>& boundary_polylines);
	static void polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points);

	Polygon_2 outer_boundary;
	std::vector<Polygon_2> holes;
};