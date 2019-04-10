#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_polyline_traits_2.h>

//#include <map.h>
using namespace std;

OrthogonalGrid::OrthogonalGrid(const CGAL::Bbox_2& bbox, int x_res, int y_res) : x_res(x_res), y_res(y_res), bbox(bbox) {
	// Create x line and y coodrinates
	create_spaced_range(bbox.xmin(), bbox.xmax(), x_res, x_coords);
	create_spaced_range(bbox.ymin(), bbox.ymax(), y_res, y_coords);
}

void OrthogonalGrid::add_additional_grid_points(const std::vector<Point_2>& additional_grid_points) {
  for (auto pt: additional_grid_points) {
    if (is_point_on_grid(pt)) continue;
    added_x_coords.push_back(pt.x());
    added_y_coords.push_back(pt.y());  
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
  std::sort(added_coords.begin(), added_coords.end());
  int start_i = 1;
  for (auto added_pt : added_coords) {
      // find the index of the closest point (but don't replace the first or last point)
      double closest_dist = std::numeric_limits<double>::infinity(); int vertex_i = -1;
      for (int i = start_i; i < axis_vec.size()-1; i++) {
        int dist = std::abs(CGAL::to_double(axis_vec[i] - added_pt));
        if (dist < closest_dist) {closest_dist = dist; vertex_i = i;}
        else break;
      }
      //std::cout << "Closest point to " << CGAL::to_double(added_pt) << " is " << axis_vec[vertex_i] << " with index = " << vertex_i << std::endl;
      // set this index as the vertex
      axis_vec[vertex_i] = added_pt;
      // now fix backwards all of the values
      Number_type min_pt = axis_vec[start_i-1];
      Number_type spacing = (added_pt-min_pt)/(vertex_i-start_i+1);
      //std::cout << "Using spacing = " << spacing << std::endl;
      //std::cout << "min_pt = " << min_pt << std::endl;
      int spacing_idx = 1;
      for (int i = start_i; i < vertex_i; i++) {
        //std::cout << "Setting index i =  " << i << " with " <<  CGAL::to_double(min_pt+spacing_idx*spacing) << " instead of " << axis_vec[i] << std::endl;
        axis_vec[i] = min_pt+spacing_idx*spacing;
        spacing_idx++;
      }
      // start the next search with another index
      start_i = vertex_i+1;
      //std::cout << "next start_i = " << start_i << std::endl;
  }
  
  // Make the rest regular
  Number_type min_pt = axis_vec[start_i-1]; int last_idx = axis_vec.size()-1; Number_type last_pt = axis_vec[last_idx];
  Number_type spacing = (last_pt-min_pt)/(last_idx-start_i+1);
  //std::cout << "min_pt = " << min_pt << std::endl;
  int spacing_idx = 1;
  for (int i = start_i; i < last_idx; i++) {
    //std::cout << "Setting index i =  " << i << " with " <<  CGAL::to_double(min_pt+spacing_idx*spacing) << " instead of " << axis_vec[i] << std::endl;
    axis_vec[i] = min_pt+spacing_idx*spacing;
    spacing_idx++;
  }
  //std::cout << "last pt = " << last_pt << std::endl;
}
void OrthogonalGrid::initialize_grid() {
  regularize_grid();
  // Add x and y lines to arrangements
  std::vector<Segment_2> grid_segments;
  for (auto x : x_coords) {grid_segments.push_back(Segment_2(Point_2(x,bbox.ymin()),Point_2(x,bbox.ymax())));}
  for (auto y : y_coords) {grid_segments.push_back(Segment_2(Point_2(bbox.xmin(),y),Point_2(bbox.xmax(),y)));}
  add_segments(grid_segments);
/*
  for (auto x : x_coords) std::cout << "x = " << x <<std::endl;
  for (auto y : y_coords) std::cout << "y = " << y <<std::endl;
  exit(1);
*/
}

Polyline_2 OrthogonalGrid::single_polyline_to_segments_on_grid(const Polyline_2& polyline, bool closed_poly) {
  std::cout << "snapping polyline to grid with closed_poly = " << closed_poly << std::endl;
   std::vector<Segment_2> grid_x_segments, grid_y_segments;
  for (auto x : x_coords) {grid_x_segments.push_back(Segment_2(Point_2(x,bbox.ymin()),Point_2(x,bbox.ymax())));}
  for (auto y : y_coords) {grid_y_segments.push_back(Segment_2(Point_2(bbox.xmin(),y),Point_2(bbox.xmax(),y)));}
  Geom_traits_2 geom_traits_2;
  std::vector<Point_2> new_poly_points;
  // Get the first polyline point
  Point_2 pt(polyline.subcurves_begin()->source());
  if (is_point_on_grid(pt)) new_poly_points.push_back(pt);
  else if (!closed_poly) {cout << "Error in OrthogonalGrid::single_polyline_to_segments_on_grid- first point " << pt << " not on grid" << endl; exit(1);}
  // If the polygon is closed we don't need the first point on the grid
  int x_coord,y_coord; x_coord = get_coord_range(pt.x(), x_coords); y_coord = get_coord_range(pt.y(), y_coords);
  Point_2 prevPt(pt);

  // Go through all polyline points
  for (auto it = polyline.subcurves_begin(); it != polyline.subcurves_end(); it++) {
    pt = it->target();
    
    //std::cout << "number_of_segments = " << it->number_of_subcurves() << std::endl;
    int new_x_coord(get_coord_range(pt.x(), x_coords)),new_y_coord(get_coord_range(pt.y(), y_coords));
    
    //cout << "pt = " << pt << " new_x_coord = " << new_x_coord << " new_y_coord = " << new_y_coord << std::endl;
    // Detect intersection
    if ( (x_coord != new_x_coord) || (y_coord != new_y_coord) ) {
      std::cout << "Intersection found from " << prevPt << " to " << pt << std::endl;
      std::cout << "\t From coords (" << x_coord << "," << y_coord << ") to (" << new_x_coord << "," << new_y_coord << ")" << std::endl;

      std::vector<Point_2> intersections;
      if (x_coord != new_x_coord) {
        std::vector<Segment_2> int_segments; int_segments.push_back(Segment_2(prevPt,pt));
        int_segments.insert( int_segments.end(), grid_x_segments.begin(), grid_x_segments.end() );
        CGAL::compute_intersection_points(int_segments.begin(), int_segments.end(),
                                  std::back_inserter(intersections), false, geom_traits_2);
      }
      if (y_coord != new_y_coord) {
        std::vector<Segment_2> int_segments; int_segments.push_back(Segment_2(prevPt,pt));
        int_segments.insert( int_segments.end(), grid_y_segments.begin(), grid_y_segments.end() );
        CGAL::compute_intersection_points(int_segments.begin(), int_segments.end(),
                                  std::back_inserter(intersections), false, geom_traits_2);
      }
      std::unique(intersections.begin(),intersections.end());
      intersections.erase(std::remove(intersections.begin(), intersections.end(), prevPt), intersections.end());
      intersections.erase(std::remove(intersections.begin(), intersections.end(), pt), intersections.end());

      int int_n = intersections.size();
      // sort them by distances
      std::vector<std::pair<int,Number_type>> int_dist;
      for (int int_i = 0; int_i < intersections.size(); int_i++) {
        int_dist.push_back(std::pair<int,Number_type>(int_i,CGAL::squared_distance(prevPt, intersections[int_i])));
      }
      std::sort(int_dist.begin(), int_dist.end(), [](const std::pair<int,Number_type>& p1, const std::pair<int,Number_type>& p2){return p1.second < p2.second;});
      for (auto i_dist : int_dist) {
        if (CGAL::squared_distance(new_poly_points.back(),intersections[i_dist.first]) < Number_type(1e-16)) {
            std::cout << "here with " << new_poly_points.back() << " and " << intersections[i_dist.first] << std::endl;
            bool are_equal  = new_poly_points.back() == intersections[i_dist.first];
            std::cout << "new_poly_points.back() == intersections[i_dist.first] is " << are_equal << std::endl;
            int b; std::cin >> b;
        } else {
          new_poly_points.push_back(intersections[i_dist.first]);
        }        
      }

      std::cout << "added pts: "; for (auto pt: intersections) {std::cout << pt <<",";} std::cout <<endl;
      /*
      if (int_n == 1) {
        new_poly_points.push_back(intersections[0]);
      } else if(int_n == 2) {
        Point_2 pt1(intersections[0]), pt2(intersections[1]);
         if (CGAL::squared_distance(prevPt, intersections[0]) > CGAL::squared_distance(prevPt, intersections[1])) std::swap(pt1,pt2);
         new_poly_points.push_back(pt1); new_poly_points.push_back(pt2);
         //std::cout << "pt1 = " << pt1 << " pt2 = " << pt2 << std::endl;
      } else if (int_n != 0) {
        std::cout << "Error at single_polyline_to_segments_on_grid: Intersections number needs to be 1 or 2 but is " << intersections.size() << std::endl;
        std::cout << "Intersections:" << std::endl;
        for (auto seg_int : intersections) std::cout << seg_int << endl;
        exit(1);
      }
      */
    }
    //std::cout << "checking if " << pt << " is on the grid" << std::endl;
    // Checking that the last point added is not the same as well (can happen when snapping singularities)
    if ( (is_point_on_grid(pt)) && (new_poly_points.back()!= pt) ) {
      new_poly_points.push_back(pt);
      std::cout << "adding pt on grid: " << pt << std::endl;
    }
    x_coord = new_x_coord; y_coord = new_y_coord;
    prevPt = pt;
  }
  int wait; std::cin >> wait;
  if ( closed_poly && (new_poly_points[0]!=new_poly_points.back()) ) new_poly_points.push_back(new_poly_points[0]);
  //cout<<"Points:"<<endl;for (auto pt: new_poly_points) std::cout << pt << endl; int wait; std::cin >> wait;
  //Polygon_2 poly2(new_poly_points.begin(),new_poly_points.end()-1); std::cout << "poly.is_simple() = " << poly2.is_simple()<< std::endl;
  //std::cout << "poly_area = " << poly2.area() << std::endl;
  //std::cout << "Ran the new one" << std::endl; int wait; cin >> wait;
  Geom_traits_2 traits;
  Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
  return polyline_construct(new_poly_points.begin(), new_poly_points.end());
}

Polyline_2 OrthogonalGrid::single_polyline_to_segments_on_grid_old(const Polyline_2& polyline) {
  //single_polyline_to_segments_on_grid2(polyline);
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
	if (!is_on_v) {std::cout << "Error! Point not on grid!" << std::endl; exit(1);}
  std::cout << "New poly" << std::endl << std::endl;
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
    // CGAL guy said they might add this functionality, but for now we need to arrange it
    // The following code just does that by tracing the curve along the new arrangement

    // We first the input curve's induces segments and save the points degree 
  	std::map<Point_2,int> point_to_deg; std::map<std::pair<Point_2,Point_2>,bool> seg_on_curve;
  	for ( auto input_edge = grid_arr.get_arrangement_internal()->induced_edges_begin(ocit); 
  			input_edge != grid_arr.get_arrangement_internal()->induced_edges_end(ocit); input_edge++) {
      auto curve_first_pt = (*input_edge)->curve();//subcurves_begin()->source();
      /*
      std::vector<Segment_2> edge_polyline_segments((*input_edge)->curve().subcurves_begin(),(*input_edge)->curve().subcurves_end());

      for (auto edge_seg: edge_polyline_segments) {
        //auto src = (*input_edge)->source(), target = (*input_edge)->target();
        //auto seg = Segment_2(src->point(),target->point());

        std::cout << "Edge = " << edge_seg << std::endl;
        auto src = edge_seg.source(); auto target = edge_seg.target();

        point_to_deg[src.point()] = src->degree();
        point_to_deg[target->point()] = target->degree();
        seg_on_curve[std::pair<Point_2,Point_2>(src->point(),target->point())] = true;
        seg_on_curve[std::pair<Point_2,Point_2>(target->point(),src->point())] = true;
      }
      */
      std::vector<Segment_2> edge_polyline_segments((*input_edge)->curve().subcurves_begin(),(*input_edge)->curve().subcurves_end());
      auto src = (*input_edge)->source(), target = (*input_edge)->target();
      auto seg = Segment_2(src->point(),target->point());

      std::cout << "Edge = " << seg << std::endl;
      std::cout << "edge_polyline_segments.size() = " << edge_polyline_segments.size() << std::endl;

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
            std::cout << "looking at point " << u->point() << std::endl;
            std::cout << "point_to_deg[u->point()]  = " << point_to_deg[u->point()]  << std::endl;
            if (point_to_deg[u->point()] > 2) {
              std::cout << "adding point " << u->point() << std::endl;
              new_poly_points.push_back(u->point());
            }
            
        }
      } while ((++curr != first) && (!found_next_pt));
      

    } while (found_next_pt);
    if (p1.x() == Number_type(301.5)) {
      std::cout << "p1  = " << p1 << " p2 = " << p2 << std::endl; int wait; std::cin >> wait;
    }
    cout<<"Points:"<<endl;for (auto pt: new_poly_points) std::cout << pt << endl;
    std::cout << "Ran the old one" << std::endl; int wait; cin >> wait;
    Geom_traits_2 traits;
    Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
    return polyline_construct(new_poly_points.begin(), new_poly_points.end());
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
  //auto it = std::lower_bound(vec.begin(),vec.end(),pt);
  //return (it != vec.end() && *it == pt);
  // works up to epsilon (rounding error from CGAL to double)
  double eps = 1e-10;
  for (auto val : vec) {
     if (std::abs(CGAL::to_double(pt-val)) < eps) return true;
  }
  return false;
}

int OrthogonalGrid::get_coord_range(const Number_type& pt, const std::vector<Number_type>& coords) {
  if (pt < coords[0]) return 0;
  if (pt > coords.back()) return coords.size()+1;
  // The point is somewhere along the grid, find the edge indices
  auto upper_bound = std::upper_bound(coords.begin(),coords.end(),pt);
  if (upper_bound == coords.end()) return coords.size();
  return upper_bound-coords.begin();
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