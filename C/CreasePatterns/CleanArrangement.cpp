#include "CleanArrangement.h"

void snap_to_boundary(const PatternBoundary& patternBoundary, const CGAL::Bbox_2& bbox, 
			PlanarArrangement& inArrangement, PlanarArrangement& outArrangement) {
	// Get all segments
	std::vector<Segment_2> inSegments; inArrangement.get_segments(inSegments);
	// Create a new list of segments without points outside the boundary, and with snapped vertices
	std::vector<Segment_2> outSegments;

	// First go through all of the vertices, and snap them to the boundary
	Arrangement_2* inArrInt = inArrangement.get_arrangement_internal();

	Arrangement_2::Vertex_const_iterator vit;
	std::cout << inArrInt->number_of_vertices() << " vertices:" << std::endl;
  	for (vit = inArrInt->vertices_begin(); vit != inArrInt->vertices_end(); ++vit) {
    	std::cout << "(" << vit->point() << ")";
    	if ((!vit->is_isolated())) {
    		// Check if it is on the boundary, and if not do something about it
    		if (!patternBoundary.on_boundary(vit->point())) {
    			// snap it to boundary, and save the snapped result in a dictionary
    		}
    	}
	}

	// Build outSegments by going through every inSegment
	// If any point is outside the polygon, don't add that segment. 
	// Otherwise, add a segment. Snap the segment vertex if necessary by the previous snapped dictionary.

	outArrangement.add_segments(outSegments);
}