
#pragma once

#include <Eigen/Dense>

#include "ArrangementDefs.h"
#include "PlanarArrangement.h"

class PatternBoundary {
public:
	PatternBoundary(std::vector<Polygon_2> boundary_polygons);
	void cut_and_snap_polylines_to_boundary(std::vector<Polyline_2>& polylines);
private:
	Polyline_2 cut_and_snap_single_polyline_to_boundary(Polyline_2& polyline);
	bool inside_boundary(const Point_2& pt);

	Polygon_2 outer_boundary;
	std::vector<Polygon_2> holes;
};