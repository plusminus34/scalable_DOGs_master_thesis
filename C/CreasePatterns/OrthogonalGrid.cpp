#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_polyline_traits_2.h>

//#include <map.h>

OrthogonalGrid::OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res) : x_res(x_res), y_res(y_res), bbox(bbox) {
	// Create x line and y coodrinates
	create_spaced_range(bbox.xmin(), bbox.xmax(), x_res, x_coords);
	create_spaced_range(bbox.ymin(), bbox.ymax(), y_res, y_coords);
}

void OrthogonalGrid::add_additional_grid_points(const std::vector<Point_2>& additional_grid_points) {
  for (auto pt : additional_grid_points) {
    subdivide_grid_at_pt(pt);
    //added_vertices.push_back(pt);
  }
}

void OrthogonalGrid::regularize_grid() {
  // Do something only if there are vertices
  if (added_x_coords.size()) {
    add_vertices_to_axis_and_keep_as_regular_as_possible(x_coords, added_x_coords);
    added_x_coords.clear();
  }
  if (added_y_coords.size()) {
    add_vertices_to_axis_and_keep_as_regular_as_possible(y_coords, added_y_coords);
    added_y_coords.clear();
  }
  //exit(1);
}
void OrthogonalGrid::add_vertices_to_axis_and_keep_as_regular_as_possible(std::vector<Number_type>& axis_vec,
  std::vector<Number_type>& added_coords) {

  int start_i = 1;
  for (auto added_pt : added_coords) {
      // find the index of the closest point (but don't replace the first or last point)
      double closest_dist = std::numeric_limits<double>::infinity(); int vertex_i = -1;
      for (int i = start_i; i < axis_vec.size()-1; i++) {
        int dist = std::abs(CGAL::to_double(axis_vec[i] - added_pt));
        if (dist < closest_dist) {closest_dist = dist; vertex_i = i;}
        else break;
      }
      std::cout << "Closest point to " << CGAL::to_double(added_pt) << " is " << axis_vec[vertex_i] << " with index = " << vertex_i << std::endl;
      // set this index as the vertex
      axis_vec[vertex_i] = added_pt;
      // now fix backwards all of the values
      Number_type min_pt = axis_vec[start_i-1];
      Number_type spacing = (added_pt-min_pt)/(vertex_i-start_i+1);
      std::cout << "min_pt = " << min_pt << std::endl;
      int spacing_idx = 1;
      for (int i = start_i; i < vertex_i; i++) {
        std::cout << "Setting index i =  " << i << " with " <<  CGAL::to_double(min_pt+spacing_idx*spacing) << " instead of " << axis_vec[i] << std::endl;
        axis_vec[i] = min_pt+spacing_idx*spacing;
        spacing_idx++;
      }
      // start the next search with another index
      start_i = vertex_i+1;
      std::cout << "next start_i = " << start_i << std::endl;
  }
  
  // Make the rest regular
  Number_type min_pt = axis_vec[start_i-1]; int last_idx = axis_vec.size()-1; Number_type last_pt = axis_vec[last_idx];
  Number_type spacing = (last_pt-min_pt)/(last_idx-start_i+1);
  std::cout << "min_pt = " << min_pt << std::endl;
  int spacing_idx = 1;
  for (int i = start_i; i < last_idx; i++) {
    std::cout << "Setting index i =  " << i << " with " <<  CGAL::to_double(min_pt+spacing_idx*spacing) << " instead of " << axis_vec[i] << std::endl;
    axis_vec[i] = min_pt+spacing_idx*spacing;
    spacing_idx++;
  }
  std::cout << "last pt = " << last_pt << std::endl;
  //exit(1);
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

	// Find the edge of the curve by going through all edges emenating from that vertex and checking for collinearity with the input
	Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
  	first = curr = v->incident_halfedges();
  	//std::cout << "Finding the edge of the polyline starting at (" << v->point() << "):";
  	Arrangement_2::Halfedge_const_handle poly_edge;
  	do {
    	// Note that the current halfedge is directed from u to v:
    	//Arrangement_2::Vertex_const_handle u = curr->source();
    	Arrangement_2::Halfedge_const_handle edge_handle = curr->twin();//->handle();

      auto curve_first_pt = edge_handle->curve().subcurves_begin()->source();
      int last_seg_i = edge_handle->curve().subcurves_end()-edge_handle->curve().subcurves_begin();
      auto first_seg = *(edge_handle->curve().subcurves_begin());
      Point_2 next_curve_point = first_seg.target();
      if (p1!= first_seg.source()) {
        next_curve_point = first_seg.source();
      }
      if (CGAL::collinear(p1, p2, next_curve_point)) {
        poly_edge = edge_handle;
        //std::cout << "found poly edge to = " << (*poly_edge)->source() << std::endl;
    		//std::cout << "Found the original polyline edge with source = " << edge_handle->source()->point() << " target = " << edge_handle->target()->point() << std::endl;
        break;
    	}
      
  	} while (++curr != first);
  	Arrangement_2::Originating_curve_iterator ocit = grid_arr.get_arrangement_internal()->originating_curves_begin(poly_edge);
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

      //std::cout << "Edge = " << seg << std::endl;

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
	//if (1) {
  /*
		x_coords.push_back(pt.x());
		std::sort(x_coords.begin(),x_coords.end());
	//} else {
		y_coords.push_back(pt.y());
		std::sort(y_coords.begin(),y_coords.end());
    */

    added_x_coords.push_back(pt.x());
    added_y_coords.push_back(pt.y());
	//}	
}

bool OrthogonalGrid::is_point_on_grid(const Point_2& pt) const {
    // First check that the pt is in the grid range
    if ((pt.x() < x_coords[0]) || (pt.x() > x_coords[x_coords.size()-1])) return false;
    if ((pt.y() < y_coords[0]) || (pt.y() > y_coords[y_coords.size()-1])) return false;

    bool is_on_x = pt_in_vec(x_coords,pt.x());
    bool is_on_y = pt_in_vec(y_coords,pt.y());
    bool is_on_added_x = pt_in_vec(added_x_coords, pt.x());
    bool is_on_added_y = pt_in_vec(added_y_coords, pt.y());

    return is_on_x || is_on_y || is_on_added_x || is_on_added_y;
}

bool OrthogonalGrid::pt_in_vec(const std::vector<Number_type>& vec, const Number_type& pt) const {
  auto it = std::lower_bound(vec.begin(),vec.end(),pt);
  return (it != vec.end() && *it == pt);
}

bool OrthogonalGrid::get_pt_edge_coordinates(const Point_2& pt, std::pair<Point_2,Point_2>& edge_pts, Number_type& t) const {
  if (!is_point_on_grid(pt)) return false;

  auto x_lower_bound = std::lower_bound(x_coords.begin(),x_coords.end(),pt.x());
  auto x_upper_bound = std::upper_bound(x_coords.begin(),x_coords.end(),pt.x());
  auto y_lower_bound = std::lower_bound(y_coords.begin(),y_coords.end(),pt.y());
  auto y_upper_bound = std::upper_bound(y_coords.begin(),y_coords.end(),pt.y());

  Point_2 v1,v2;
  // The point is on the grid, so either *x_lower_bound == pt.x() or *y_lower_bound == pt.y() or both (in case it's a vertex)
  if ( (*x_lower_bound == pt.x()) && (*y_lower_bound == pt.y()) ) {
    t = 1.; v1 = v2 = pt;
  } else if (*x_lower_bound == pt.x()) {
    y_lower_bound--; // get one element before (so it will be smaller)
    v1 = Point_2(pt.x(),*y_lower_bound); v2 = Point_2(pt.x(),*y_upper_bound);
    t = (pt.y()-v2.y())/(v1.y()-v2.y());
    // compute t
  } else {
    x_lower_bound--; // get one element before (so it will be smaller)
    v1 = Point_2(*x_lower_bound, pt.y()); v2 = Point_2(*x_upper_bound, pt.y());
    t = (pt.x()-v2.x())/(v1.x()-v2.x());
  }

  edge_pts = std::pair<Point_2,Point_2>(v1,v2);
  return true;
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