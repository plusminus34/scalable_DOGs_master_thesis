/*
#pragma once

#include <Eigen/Dense>

#include "ArrangementDefs.h"
#include "PlanarArrangement.h"

class PatternBoundary {
public:
	PatternBoundary(Polyline_2& polyline);
	void cut_and_snap_polylines_to_boundary(std::vector<Polyline_2>& polylines);
	Polyline_2 get_polyline() {return bnd_poly;};
private:
	Polyline_2 cut_and_snap_single_polyline_to_boundary(Polyline_2& polyline);

	Polyline_2& bnd_poly;
	PlanarArrangement arrangement;
};
*/