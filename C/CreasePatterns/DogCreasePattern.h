#pragma once

#include <Eigen/Dense>

#include "PlanarArrangement.h"
#include "OrthogonalGrid.h"


class DogCreasePattern {
  
public:
	DogCreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, int x_res, int y_res, bool snap_rounding = false);
	// TODO constructor from polygon (requiring to support removal of faces from the grid)
/*
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
*/

private:
	// snap rounding (and possibly later project initial curves to boundary)
	void init_initial_arrangement_and_polylines(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, bool snap_rounding);
	void bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline);

	std::vector<Polyline_2> initial_polylines; // The boundary is the first one
	PlanarArrangement initial_arrangement;
	OrthogonalGrid orthogonalGrid;
	//PlanarArrangement clipped_grid_arrangement;
};