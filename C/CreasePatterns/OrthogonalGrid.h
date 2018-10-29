#pragma once

#include <Eigen/Dense>

#include "PlanarArrangement.h"


class OrthogonalGrid {
  
public:
	// bounding box and resolution. Set up arrangement (have one as member)
	OrthogonalGrid(){};

	// multiple polylines
	void polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines);
private:
	PlanarArrangement gridArrangement;
};
