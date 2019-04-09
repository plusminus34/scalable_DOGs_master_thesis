#pragma once

#include <Eigen/Dense>

#include <CGAL/Bbox_2.h>
#include "ArrangementDefs.h"

#include "PlanarArrangement.h"

class OrthogonalGrid : public PlanarArrangement {
  
public:
	// x_res are number of x vertices, y_res number of y vertices (so number_of_edges-1)
	// bounding box and resolution. Set up arrangement (have one as member)
	// Third parameter is used for curves intersection, hence we need to subdivide the grid to add that
	OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res);
	void add_additional_grid_points(const std::vector<Point_2>& additional_grid_points);

	void initialize_grid();

	// multiple polylines
	//std::vector<Polyline_2> polylines_to_segments_on_grid(const std::vector<Polyline_2>& polylines);
	Polyline_2 single_polyline_to_segments_on_grid(const Polyline_2& polyline, bool closed_poly = false);
	Polyline_2 single_polyline_to_segments_on_grid_old(const Polyline_2& polyline);

	const std::vector<Number_type>& get_x_coords() const {return x_coords;}
	const std::vector<Number_type>& get_y_coords() const {return y_coords;}

	// Return true if the point pt is on an edge, false otherwise
	// In case it is on an edge, return the two vertices next to it v1,v2, and a weight t s.t pt = t*v1 + (1-t)*v2
	// In case it is on a vertex, t will be 1 (and v1=v2 but in that case one should ignore v2)
	bool get_pt_edge_coordinates(const Point_2& pt, std::pair<Point_2,Point_2>& edge_pts, Number_type& t) const;

	Number_type get_avg_x_edge() const {return (bbox.xmax()-bbox.xmin())/(x_res-1);}
	Number_type get_avg_y_edge() const {return (bbox.ymax()-bbox.ymin())/(y_res-1);}

private:
	void regularize_grid();
	// True if it is on the grid lines (edges or vertices)
	void add_vertices_to_axis_and_keep_as_regular_as_possible(std::vector<Number_type>& axis_vec, std::vector<Number_type>& added_coords);
	bool is_point_on_grid(const Point_2& pt) const;
	void create_spaced_range(const Number_type min, const Number_type max, const int num_points, std::vector<Number_type>& range);

	void polyline_to_segments(const Polyline_2& polyline, std::vector<Segment_2>& segments);
	bool pt_in_vec(const std::vector<Number_type>& vec, const Number_type& pt) const;

	static int get_coord_range(const Number_type& pt, const std::vector<Number_type>& coords);
	//Number_type lower_bound(const Point_2& pt, int axis);
	// Holds only the grid, as there's no need to have multiple intersected polylines calculation all the time
	
	double x_res,y_res;
	const CGAL::Bbox_2& bbox;
	std::vector<Number_type> added_x_coords, added_y_coords;
	std::vector<Number_type> x_coords;
	std::vector<Number_type> y_coords;
};
