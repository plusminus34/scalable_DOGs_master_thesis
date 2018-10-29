#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


void OrthogonalGrid::polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines) {
	// Compute the non-intersecting sub-segments induced by the input segments.
	//std::vector<Polyline_2> sub_polylines;
	//CGAL::compute_subcurves(polylines.begin(), polylines.end(), std::back_inserter(sub_polylines));
	std::vector<Segment_2> in, out;
	CGAL::compute_subcurves(polylines.begin(), polylines.end(), std::back_inserter(out));
	

}