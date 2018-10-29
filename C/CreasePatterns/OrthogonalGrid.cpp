#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


void OrthogonalGrid::polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines) {
	// Compute the non-intersecting sub-segments induced by the input segments.
	Geom_traits_2 geom_traits_2;
	std::vector<Polyline_2_Monotone> sub_polylines;
	CGAL::compute_subcurves(polylines.begin(), polylines.end(), std::back_inserter(sub_polylines), false, geom_traits_2);

	// We expect here to have the sub_polylines start and end segment touching the grid. In case they don't, this is a problem
}