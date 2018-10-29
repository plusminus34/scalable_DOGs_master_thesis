#include "PatternBoundary.h"


PatternBoundary::PatternBoundary(Polyline_2& polyline) : bnd_poly(polyline) {
	arrangement.add_polyline(polyline);
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