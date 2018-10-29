#include "PlanarArrangement.h"

void PlanarArrangement::add_polylines(std::vector<Polyline_2>& polylines) {
	insert(arr, polylines.begin(), polylines.end());
}

void PlanarArrangement::add_polyline(Polyline_2& polyline) {
	insert(arr, polyline);
}