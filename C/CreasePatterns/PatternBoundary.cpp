
#include "PatternBoundary.h"


PatternBoundary::PatternBoundary(const std::vector<Polygon_2>& boundary_polygons) {
	std::vector<Number_type> areas; Number_type max_area(0); int boundary_poly_i = 0;
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (!boundary_polygons[i].is_simple()) {
			std::cout << "Error at PatternBoundary: Polygon " << i << " is not simple" << std::endl;
			exit(1);
		}
		if (abs(boundary_polygons[i].area()) > max_area) {
			max_area = abs(boundary_polygons[i].area());
			boundary_poly_i = i;
		}
	}
	outer_boundary = boundary_polygons[boundary_poly_i];
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (i!= boundary_poly_i) {holes.push_back(Polygon_2(boundary_polygons[i]));}
	}
}

PatternBoundary::PatternBoundary(const std::vector<Polyline_2>& boundary_polylines) : 
	PatternBoundary(polylines_to_polygons(boundary_polylines)) {
	// empty on purpose
}

Polyline_2 PatternBoundary::clip_and_snap(Polyline_2& polyline) {
	// Get polyline points
	return polyline; // todo implement
}

bool PatternBoundary::inside_boundary(const Point_2& pt) {
	// True if the point is inside the outer boundary but outside the holes (and including the boundary of both)
	if (outer_boundary.bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE) return false;
	for (int i = 0; i < holes.size(); i++) if (holes[i].bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) return false;
	return true;
}


std::vector<Polygon_2> PatternBoundary::polylines_to_polygons(const std::vector<Polyline_2>& boundary_polylines) {
	std::vector<Polygon_2> polygons;
	for (auto poly:boundary_polylines) {
		std::vector<Point_2> points; polyline_to_points(poly, points);
		// Make sure the first and last point are not equal
		if (points[0] == points.back()) points.pop_back();
		polygons.push_back(Polygon_2(points.begin(),points.end()));
	}
	return polygons;
}

void PatternBoundary::polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points) {
	int seg_n = poly.subcurves_end()-poly.subcurves_begin();
	points.push_back(poly.subcurves_begin()->source());
	for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
		points.push_back(seg_i->target());
	}
}