#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_polyline_traits_2.h>

//#include <map.h>

OrthogonalGrid::OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res, const std::vector<Point_2>& additional_grid_points) {
	// Create x line and y coodrinates
	create_spaced_range(bbox.xmin(), bbox.xmax(), x_res, x_coords);
	create_spaced_range(bbox.ymin(), bbox.ymax(), y_res, y_coords);
	for (auto pt : additional_grid_points) {subdivide_grid_at_pt(pt);}
	// Add x and y lines to arrangements
	std::vector<Segment_2> grid_segments;
	for (auto x : x_coords) {grid_segments.push_back(Segment_2(Point_2(x,bbox.ymin()),Point_2(x,bbox.ymax())));}
	for (auto y : y_coords) {grid_segments.push_back(Segment_2(Point_2(bbox.xmin(),y),Point_2(bbox.xmax(),y)));}
	add_segments(grid_segments);
}

void OrthogonalGrid::polylines_to_segments_on_grid(std::vector<Polyline_2>& polylines) {
	for (auto poly: polylines) {
		poly = single_polyline_to_segments_on_grid(poly);
	}
}


// Todo add polylines
// Todo snap polyline to original one
//		Find the initial vertex on the polyline and from there traverse the original polyline (find the curve by next vertex that should be first segment..)
//		then check the degree of every vertex
Polyline_2 OrthogonalGrid::single_polyline_to_segments_on_grid(const Polyline_2& polyline) {
	// Create a copy of the grid arrangement and add the polyline to it
	PlanarArrangement grid_arr(*this);
	grid_arr.add_polyline(polyline);
	std::vector<Segment_2> poly_initial_segments; polyline_to_segments(polyline, poly_initial_segments);
	std::vector<Point_2> new_poly_points;

	// Insert the initial vertex (assume it is on the grid!)
	Point_2 p1(poly_initial_segments[0].source()), p2(poly_initial_segments[0].target());
	// Find the initial vertex on the polyline
	Vertex_const_handle v;
	bool is_on_v = grid_arr.locate_point_on_vertex(p1,v);
	if (!is_on_v) std::cout << "Error! Point not on grid!" << std::endl;

	// Go through all edges emenating from v
	Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
  	first = curr = v->incident_halfedges();
  	std::cout << "Finding the edge of the polyline starting at (" << v->point() << "):";
  	Arrangement_2::Halfedge_const_handle poly_edge;
  	do {
    	// Note that the current halfedge is directed from u to v:
    	Arrangement_2::Vertex_const_handle u = curr->source();
    	Arrangement_2::Halfedge_const_handle edge_handle = curr->twin();//->handle();
    	// Todo: In case the original edge is and edge of the grid, this could fail, and we will need further checks (origin of the edge)
    	//	This is possible with the history data but not needed atm
    	if (CGAL::collinear(p1, p2, edge_handle->target()->point())) {
    		poly_edge = edge_handle;
    		std::cout << "Found the original polyline edge with source = " << edge_handle->source()->point() << " target = " << edge_handle->target()->point() << std::endl;
    	}
  	} while (++curr != first);
  	Arrangement_2::Originating_curve_iterator ocit = grid_arr.get_arrangement_internal()->originating_curves_begin(poly_edge);
  	// For now assume that this edge as only one originating curve (todo this can be false)
  	// TODO in the other case (degenerate polyline on grid), iterate with
  	/*
  	for (ocit = grid_arr.get_arrangement_internal()->originating_curves_begin(eit); 
  				ocit != grid_arr.get_arrangement_internal()->.originating_curves_end(eit); ++ocit) {
  	*/
  	// 
  	std::cout << "originating curve is " << *ocit << std::endl;
  	//Arrangement_2::Halfedge_handle input_edge;
  	//input_edge = grid_arr.get_arrangement_internal()->induced_edges_begin(ocit);
  	
  	/*
  	for ( auto input_edge = grid_arr.get_arrangement_internal()->induced_edges_begin(ocit); input_edge != grid_arr.get_arrangement_internal()->induced_edges_end(ocit); input_edge++) {
  		//std::cout << (*input_edge)->source()->point() << std::endl;
  		std::cout << "Found induced edge with source = " << (*input_edge)->source()->point() << " and target = " << (*input_edge)->target()->point() << std::endl;
  		std::cout << "The target degree is " << (*input_edge)->target()->degree() << std::endl;
  		if ((*input_edge)->target()->degree() > 2) {
  			new_poly_points.push_back((*input_edge)->target()->point());
  		}
  	}
  	*/
  	std::vector<std::pair<int,int>> seg_deg;
  	std::vector<Segment_2> unsorted_segments; //std::vector<Arrangement_2::Halfedge_handle> edge_handle_vec;
  	for ( auto input_edge = grid_arr.get_arrangement_internal()->induced_edges_begin(ocit); 
  			input_edge != grid_arr.get_arrangement_internal()->induced_edges_end(ocit); input_edge++) {
		unsorted_segments.push_back(Segment_2((*input_edge)->source()->point(),(*input_edge)->target()->point()));
		//std::cout << "degree = " << (*input_edge)->target()->degree() << std::endl;
		seg_deg.push_back(std::pair<int,int>((*input_edge)->source()->degree(),(*input_edge)->target()->degree()));
		//edge_handle_vec.push_back(*input_edge);
  	}
  	// TODO: The following is kinda disgusting. There gotta be a better way to do that but this is due the edges of induced_edges_begin() 
  	//			not being sorted by the curve. Otherwise I'll write a better faster manual algorithm (single_polyline_to_segments_on_grid_fast) 
  	//			because this makes no sense.
  	// This is still unsorted, just make the search much more efficient
  	// std::sort(unsorted_segments.begin(),unsorted_segments.end());
  	// Sadly the segments can even be reversed
  	Segment_2 last_seg(poly_edge->source()->point(), poly_edge->target()->point());
  	new_poly_points.push_back(last_seg.source());//new_poly_points.push_back(last_seg.target());
  	int last_point_idx = 0;
  	for (int i = 0; i < unsorted_segments.size()-1; i++) {
  		// Find a segment with source as polyline_segments[polyline_segments]
  		Point_2 next_point = last_seg.target();
  		
  		bool add_point = false;
  		//for (auto seg : unsorted_segments) {
  		for (int j = 0; j < unsorted_segments.size(); j++) {
  			auto seg = unsorted_segments[j];
  			// Possible continuation of the segment
  			if ((last_seg.target() == seg.source()) || (last_seg.target() == seg.target()) ) {
  				// Not the same segment
  				if ((seg.source()!= last_seg.source()) || (seg.target()!= last_seg.target()))  {
  					// If the next segment is a continuation but a flip, just flip it
  					if (next_point == seg.target()) {
  						seg = Segment_2(seg.target(),seg.source());
  						add_point = (seg_deg[j].first > 2);
  					} else {
  						add_point = (seg_deg[j].second > 2);
  					}
  					last_seg = seg;
  					//add_point = true;
  					if (add_point) {
  						new_poly_points.push_back(seg.target());
  						last_point_idx++;
  					}
  					continue;
  				}
  			}
  		}
  	}

	// Construct the polyline from the points
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	Polyline_2 pi = polyline_construct(new_poly_points.begin(), new_poly_points.end());
	return pi;
}

void OrthogonalGrid::subdivide_grid_at_pt(const Point_2& pt) {
	if (is_point_on_grid(pt)) return;
	// Todo: check what is better x or y
	if (1) {
		x_coords.push_back(pt.x());
		std::sort(x_coords.begin(),x_coords.end());
	} else {
		x_coords.push_back(pt.y());
		std::sort(y_coords.begin(),y_coords.end());
	}	
}

bool OrthogonalGrid::is_point_on_grid(const Point_2& pt) {
	auto it = std::lower_bound(x_coords.begin(),x_coords.end(),pt.x(),
        [](const Number_type& l, const Number_type& r){ return l < r; });
    bool is_on_x = (it != x_coords.end() && *it == pt.x());

    it = std::lower_bound(y_coords.begin(),y_coords.end(),pt.y(),
        [](const Number_type& l, const Number_type& r){ return l < r; });
    bool is_on_y = (it != y_coords.end() && *it == pt.y());

    return is_on_x || is_on_y;
}

void OrthogonalGrid::create_spaced_range(const Number_type min, const Number_type max, const int num_points, std::vector<Number_type>& range) {
	Number_type spacing = (max-min)/(num_points-1);
	for (int i = 0; i < num_points; i++) {range.push_back(min+i*spacing);}
}

void OrthogonalGrid::polyline_to_segments(const Polyline_2& polyline, std::vector<Segment_2>& segments) {
	for (auto it = polyline.subcurves_begin(); it != polyline.subcurves_end(); it++) {
		segments.push_back(*it);
	}
}

/*
Polyline_2 OrthogonalGrid::single_polyline_to_segments_on_grid_fast(const Polyline_2& polyline) {
	// NOT IMPLEMENTED (and could be used probably just for rectangular grids boundary without additional changes)
	std::vector<Point_2> grid_poly;
	std::vector<Segment_2> segments; polyline_to_segments(polyline, segments);

	// Assume first points lies on edges
	Point_2 pt = segments[0].source(); int x_i, y_i;
	grid_poly.push_back(pt);


	// 1) Add all intersections of this edge. Should recive a given face (the one we are in at the beginning). 
	//		We need to start by testing all edges around the face besides the one we already know we intersect
	//		Then we enter a new face index, and need again to consider all intersecting edges in the face, besides the one we came from
	//		We can just filter the previous intersection by not adding grid_poly[-1] in python notation (the last number)
	//		So we need to know which face are we in, and get 4 segments around it.
	//		The routine at the end should return the face of the last segment.
	//		We should handle the case where we lie exactly on a vertex (then we should check four faces)
	//		This should stop once we don't find any intersections.

	// Repeat the process when we update our current face at every point to check for edges.
	


	//locate_grid_xy_indices()

	// Todo: implement
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	std::list<Point_2> points;
	Polyline_2 pi = polyline_construct(points.begin(), points.end());
	return pi;
}
*/