
#include "PatternBoundary.h"


PatternBoundary::PatternBoundary(std::vector<Polygon_2> boundary_polygons) {
	std::vector<Number_type> areas; Number_type max_area(0); int boundary_poly_i = 0;
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (boundary_polygons[i].area() > max_area) {
			max_area = boundary_polygons[i].area();
			boundary_poly_i = i;
		}
	}
	outer_boundary = boundary_polygons[boundary_poly_i];
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (i!= boundary_poly_i) {holes.push_back(Polygon_2(boundary_polygons[i]));}
	}
}

void PatternBoundary::cut_and_snap_polylines_to_boundary(std::vector<Polyline_2>& polylines) {
	for (auto poly : polylines) {
		//poly = cut_and_snap_single_polyline_to_boundary(poly);
	}
}

Polyline_2 PatternBoundary::cut_and_snap_single_polyline_to_boundary(Polyline_2& polyline) {
	// Get polyline points
	return polyline; // todo implement
}

bool PatternBoundary::inside_boundary(const Point_2& pt) {
	// True if the point is inside the outer boundary but outside the holes (and including the boundary of both)
	if (outer_boundary.bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE) return false;
	for (int i = 0; i < holes.size(); i++) if (holes[i].bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) return false;
	return true;
}
