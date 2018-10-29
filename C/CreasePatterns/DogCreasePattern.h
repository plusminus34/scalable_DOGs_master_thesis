#pragma once

#include <Eigen/Dense>

#include "PatternBoundary.h"
#include "OrthogonalGrid.h"

/*
class DogCreasePattern {
  
public:
	DogCreasePattern(BoundingBox& borderBox, std::vector<Polylines>& polylines){};

	// multiple polylines
	void get_initial_arrangement();
	void place_dog_on_arrangement(int x, int y); // also subdivide stuff
	void get_arrangement_after_grid_intersection();
private:
	void clean_and_snap_to_boundary();

	void find_polylines_as_grid_intersections(); // new polylines as intersections with the grid
	void subdivide_grid_at_intersections();
	void find_curves_intersection(); // do it by vertex degree and not being on the boundary
	void snap_rounding();
	void is_on_boundary();
};
*/