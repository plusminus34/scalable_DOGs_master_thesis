#pragma once

#include <Eigen/Dense>

#include <CGAL/Bbox_2.h>
#include "ArrangementDefs.h"

#include "PlanarArrangement.h"

class OrthogonalGrid {
  
public:
	// x_res are number of x vertices, y_res number of y vertices (so number_of_edges-1)
	// bounding box and resolution. Set up arrangement (have one as member)
	// Third parameter is used for curves intersection, hence we need to subdivide the grid to add that
	OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res, const std::vector<Point_2>& additional_grid_points);

	// multiple polylines
	void polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines);
private:
	bool is_point_on_grid(const Point_2& pt);
	void subdivide_grid_at_pt(const Point_2& pt);
	void create_spaced_range(const Number_type min, const Number_type max, const int num_points, std::vector<Number_type>& range);

	Polyline_2 single_polyline_to_segments_on_grid(const Polyline_2& polyline);

	void polyline_to_segments(const Polyline_2& polyline, std::vector<Segment_2>& segments);

	

	//Number_type lower_bound(const Point_2& pt, int axis);
	// Holds only the grid, as there's no need to have multiple intersected polylines calculation all the time
	PlanarArrangement arrangement;
	std::vector<Number_type> x_coords;
	std::vector<Number_type> y_coords;
};
