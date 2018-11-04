#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_polyline_traits_2.h>

//#include <map.h>

OrthogonalGrid::OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res) : bbox(bbox) {
	// Create x line and y coodrinates
	create_spaced_range(bbox.xmin(), bbox.xmax(), x_res, x_coords);
	create_spaced_range(bbox.ymin(), bbox.ymax(), y_res, y_coords);
}

void OrthogonalGrid::add_additional_grid_points(const std::vector<Point_2>& additional_grid_points) {
  for (auto pt : additional_grid_points) {subdivide_grid_at_pt(pt);}
}

void OrthogonalGrid::initialize_grid() {
  // Add x and y lines to arrangements
  std::vector<Segment_2> grid_segments;
  for (auto x : x_coords) {grid_segments.push_back(Segment_2(Point_2(x,bbox.ymin()),Point_2(x,bbox.ymax())));}
  for (auto y : y_coords) {grid_segments.push_back(Segment_2(Point_2(bbox.xmin(),y),Point_2(bbox.xmax(),y)));}
  add_segments(grid_segments);
}

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

  std::cout << "yo" << std::endl;
	// Find the edge of the curve by going through all edges emenating from that vertex and checking for collinearity with the input
	Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
  	first = curr = v->incident_halfedges();
  	//std::cout << "Finding the edge of the polyline starting at (" << v->point() << "):";
  	Arrangement_2::Halfedge_const_handle poly_edge;
  	do {
    	// Note that the current halfedge is directed from u to v:
    	Arrangement_2::Vertex_const_handle u = curr->source();
    	Arrangement_2::Halfedge_const_handle edge_handle = curr->twin();//->handle();


      Point_2 next_curve_point = edge_handle->curve().subcurves_begin()->target();
      if (edge_handle->curve().subcurves_begin()->source() != curr->source()->point()) {
        int last_seg_i = edge_handle->curve().subcurves_end()-edge_handle->curve().subcurves_begin();
        next_curve_point = edge_handle->curve().subcurves_begin()->source();
      }
      

    	// Todo: In case the original edge is and edge of the grid, this could fail, and we will need further checks (origin of the edge)
    	//	This is possible with the history data but not needed atm
    	//if (CGAL::collinear(p1, p2, edge_handle->target()->point())) {
      if (CGAL::collinear(p1, p2, next_curve_point)) {
    		poly_edge = edge_handle;
    		//std::cout << "Found the original polyline edge with source = " << edge_handle->source()->point() << " target = " << edge_handle->target()->point() << std::endl;
    	}
  	} while (++curr != first);
    std::cout << "da" << std::endl;
  	Arrangement_2::Originating_curve_iterator ocit = grid_arr.get_arrangement_internal()->originating_curves_begin(poly_edge);
    std::cout << "bla" << std::endl;
  	// For now assume that this edge has only one originating curve (todo this can be false)
  	
    // We now have all the edges induced by the original input, the problem is that they are not ordered.
    // CGAL guy said they might add this functionality, but for now we need to arrange tit
    // The following code just does that by tracing the curve along the new arrangement

    // We first the input curve's induces segments and save the points degree 
  	std::map<Point_2,int> point_to_deg; std::map<std::pair<Point_2,Point_2>,bool> seg_on_curve;
  	for ( auto input_edge = grid_arr.get_arrangement_internal()->induced_edges_begin(ocit); 
  			input_edge != grid_arr.get_arrangement_internal()->induced_edges_end(ocit); input_edge++) {
      auto src = (*input_edge)->source(), target = (*input_edge)->target();
      auto seg = Segment_2(src->point(),target->point());

      std::cout << "Edge from " << seg << " to " << seg << std::endl;

      point_to_deg[src->point()] = src->degree();
      point_to_deg[target->point()] = target->degree();
      seg_on_curve[std::pair<Point_2,Point_2>(src->point(),target->point())] = true;
      seg_on_curve[std::pair<Point_2,Point_2>(target->point(),src->point())] = true;
  	}

    // traverse the graph, choosing the next vertex has one in the key of the map, and adding it in case its degree > 2
    //    The latter means it intersects the grid (or another Construct_curve_2)
    Point_2 prev_pt = poly_edge->source()->point(); auto cur_pt_ptr = poly_edge->target();
    point_to_deg[prev_pt] = 0;
    new_poly_points.push_back(prev_pt);
    if (point_to_deg[cur_pt_ptr->point()]>2) new_poly_points.push_back(cur_pt_ptr->point());
    
    bool found_next_pt = false;
    do {
      // go through all the neighbours of cur_pt, check that they are marked and different then prev_pt
      Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
      first = curr = cur_pt_ptr->incident_halfedges();
      
      found_next_pt = false;
      do {
        // Note that the current halfedge is directed from u to v:
        Arrangement_2::Vertex_const_handle u = curr->source();
        std::pair<Point_2,Point_2> nb_seg(cur_pt_ptr->point(),u->point());
        
        if ((seg_on_curve.count(nb_seg)) && (point_to_deg[u->point()] > 0)) {
          
            point_to_deg[cur_pt_ptr->point()] = 0;
            cur_pt_ptr = u;
            found_next_pt = true;
            if (point_to_deg[u->point()] > 2) {
              new_poly_points.push_back(u->point());
            }
            
        }
      } while ((++curr != first) && (!found_next_pt));
      

    } while (found_next_pt);
    Geom_traits_2 traits;
    Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
    return polyline_construct(new_poly_points.begin(), new_poly_points.end());
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