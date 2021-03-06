
#pragma once

#include <Eigen/Dense>

#include "ArrangementDefs.h"
#include "PlanarArrangement.h"

class PatternBoundary {
public:
	PatternBoundary(const std::vector<Polygon_2>& boundary_polygons);
	PatternBoundary(const std::vector<Polyline_2>& boundary_polylines);
	Polyline_2 filter_and_snap(Polyline_2& polyline, const Number_type& squared_dist_threshold,
			Point_2& firstPt, Point_2& lastPt, bool& snappedFirst, bool& snappedLast);

	std::vector<Point_2> get_all_boundary_points() const;
	bool inside(const Point_2& pt);
	bool on_boundary(const Point_2& pt) const;
	std::vector<Point_2> get_vertical_and_horizontal_intersections(Point_2& pt);

	std::vector<Polygon_2> get_holes() const {return holes;}
	static void polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points);
	static void proj_pt_to_segment(const Point_2& pt, const Segment_2& seg, Point_2& proj_pt, Number_type& squared_dist);
private:
	Point_2 snap_pt(const Point_2& pt, const Number_type& squared_dist_threshold, bool& has_snapped);

	static std::vector<Polygon_2> polylines_to_polygons(const std::vector<Polyline_2>& boundary_polylines);

	///bool is_hole_up_to_tolerance(Polygon_with_hole_2)

	Polygon_2 outer_boundary;
	std::vector<Polygon_2> holes;
	std::vector<Segment_2> all_polygons_edges;
};