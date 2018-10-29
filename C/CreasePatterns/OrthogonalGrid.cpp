#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


void OrthogonalGrid::polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines) {
	CGAL::compute_intersection_points(polylines.begin(), polylines.end(),
                                    std::back_inserter(pts));

}