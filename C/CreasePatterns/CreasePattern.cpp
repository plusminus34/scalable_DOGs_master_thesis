#include "CreasePattern.h"

#include "OrthogonalGrid.h"

#include <igl/combine.h>

CreasePattern::CreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2> polylines, std::vector<Polyline_2> bnd_polylines,
							int x_res, int y_res) :
											bbox(bbox), orthogonalGrid(bbox, x_res, y_res) {
	// Threshold used for snapping
	double bbox_max_len = std::max(std::abs(CGAL::to_double(bbox.xmax()-bbox.xmin())),std::abs(CGAL::to_double(bbox.ymax()-bbox.ymin())));
	Number_type dist_threshold_pow2(pow(bbox_max_len/50,2)); Number_type is_closed_threshold(pow(bbox_max_len/500,2));

	// Handle polyline intersections
	std::vector<Point_2> endpoints_intersections;
	auto merged_starting_point_polylines = snap_nearby_polylines_start_end_starting_points(polylines, endpoints_intersections);
	// At that point just make sure to snap the 'x' and 'y' coordinates of the intersections to one of them first
	initial_fold_polylines = merge_nearby_polylines_intersections(merged_starting_point_polylines);

	// Setup initial boundary (for now just a boundary box)
	initial_bnd_polylines = bnd_polylines;

	// The following is just for visualization
	initial_arrangement.add_polylines(initial_fold_polylines); initial_arrangement.add_polylines(initial_bnd_polylines);

	// get poly lines intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(initial_fold_polylines.begin(), initial_fold_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
  	polylines_intersections.insert(polylines_intersections.end(), endpoints_intersections.begin(), endpoints_intersections.end());
  	std::unique(polylines_intersections.begin(),polylines_intersections.end());
  	// add segment start ending intersection points (passing true in the above function always gives back all endpoints even if they dont intersect)
  	std::cout << "polylines_intersections.size() = " << polylines_intersections.size() << std::endl;
  	for (auto p : polylines_intersections) {std::cout << "Intersection at " << p << std::endl;}
  	
  	// Create an orthogonal grid with singularities
  	//orthogonalGrid.add_additional_grid_points(polylines_intersections);
  	std::vector<Point_2> crease_vertices = polylines_intersections;

  	// add additional interesction points at the start and end of every curve, after snapping
  	PatternBoundary origPatternBoundary(initial_bnd_polylines);
  	std::vector<Polyline_2> filtered_and_clipped_to_boundary_polylines;
  	for (auto poly = initial_fold_polylines.begin(); poly != initial_fold_polylines.end(); poly++) {
  		bool closed_polyline = is_polyline_closed_with_tolerance(*poly, is_closed_threshold);
  		if (!closed_polyline) {
  			Point_2 firstPt,lastPt; bool snappedFirst,snappedLast;
	  		filtered_and_clipped_to_boundary_polylines.push_back(
	  				origPatternBoundary.filter_and_snap(*poly,dist_threshold_pow2, firstPt, lastPt, snappedFirst, snappedLast));
	  		// Add the first and last point to the grid
	  		if (snappedFirst) { std::cout << "adding point " << firstPt << " to grid" << std::endl; crease_vertices.push_back(firstPt);}
	  		if (snappedLast)  { std::cout << "adding point " << lastPt << " to grid" << std::endl; crease_vertices.push_back(lastPt);}
  		} else {
  			filtered_and_clipped_to_boundary_polylines.push_back(*poly);
  		}
  	}
  	crease_vertices = align_crease_vertices_x_y_with_boundary(origPatternBoundary, crease_vertices, polylines_intersections.size(), dist_threshold_pow2);
  	auto preprocessed_polylines = snap_polylines_start_end_to_vertices(filtered_and_clipped_to_boundary_polylines, crease_vertices, dist_threshold_pow2);

  	orthogonalGrid.add_additional_grid_points(crease_vertices);
  	orthogonalGrid.initialize_grid();
  	
	// Clip boundary polylines to grid
	for (auto poly = initial_bnd_polylines.begin(); poly != initial_bnd_polylines.end(); poly++) {
		bool closed_polyline = true;
		clipped_bnd_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(*poly, closed_polyline));
	}
	
	// Create boundary
	patternBoundary = new PatternBoundary(clipped_bnd_polylines);
	
	// Clip fold polylines to grid, clip and snap them to the boudnary 	
	//for (auto poly = filtered_and_clipped_to_boundary_polylines.begin(); poly != filtered_and_clipped_to_boundary_polylines.end(); poly++) {
	for (auto poly : preprocessed_polylines) {
		bool closed_polyline = is_polyline_closed_with_tolerance(poly, is_closed_threshold);
		clipped_fold_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(poly,closed_polyline));
	}
	std::cout << "clipped polylines to boundary" << std::endl;

	clipped_grid_arrangement.add_polylines(clipped_fold_polylines);
	clipped_grid_arrangement.add_polylines(clipped_bnd_polylines);
	
}

CreasePattern::CreasePattern(const CreasePattern& cP) : initial_fold_polylines(cP.initial_fold_polylines), initial_bnd_polylines(cP.initial_bnd_polylines), initial_arrangement(cP.initial_arrangement),
					orthogonalGrid(cP.orthogonalGrid), clipped_fold_polylines(cP.clipped_fold_polylines), clipped_bnd_polylines(cP.clipped_bnd_polylines), clipped_grid_arrangement(cP.clipped_grid_arrangement), 
					bbox(cP.bbox) {
	// empty
}

bool CreasePattern::get_snapped_vertices_locations(const std::vector<Point_2>& polylines_intersections, Number_type threshold,
	 std::map<Point_2, Point_2>& vertices_to_snapped_vertices) {
	std::cout << "getting snapped vertices locations" << std::endl;
  	// Find nearby intersection groups and have a map between original vertices and a new coordinate
  	bool should_snap = false;
  	std::vector<std::set<int> > indices_groups; std::vector<bool> in_group(polylines_intersections.size(), false);
  	for (int i = 0; i < polylines_intersections.size(); i++) {
  		std::cout << "vertex at: " << polylines_intersections[i] << std::endl;
  		std::set<int> nearby_points; nearby_points.insert(i);
  		for (int j = i+1; j < polylines_intersections.size(); j++) {
  			auto diff = polylines_intersections[i]-polylines_intersections[j];
  			if ( (diff.x()*diff.x()+diff.y()*diff.y()) < threshold) {
  				nearby_points.insert(j);
  				should_snap = true;
  			}
  		}
  		// see if j was already added somewhere
  		bool new_set = true;
  		for (int k = 0; k < indices_groups.size(); k++) {
  			if (indices_groups[k].count(i)) {
  				new_set = false;
  				for (auto idx: nearby_points) indices_groups[k].insert(idx);
  				continue;
  			}
  		}
  		if (new_set) indices_groups.push_back(nearby_points);
  		// 
  		/*
  		if (!in_group[i]) {
  			indices_groups.push_back()
  			in_group[i] = true;
  		}
  		*/
  	}

	for (int k = 0; k < indices_groups.size(); k++) {
		// Compute the average of the points there
		Point_2 new_v(0,0);
		for (auto pt_i : indices_groups[k]) {new_v = Point_2(new_v.x() + polylines_intersections[pt_i].x(), new_v.y() + polylines_intersections[pt_i].y());}
		new_v = Point_2(new_v.x() / Number_type(indices_groups[k].size()),new_v.y() / Number_type(indices_groups[k].size()));

		std::cout << "Vertices group " << k << ": "; 
		for (auto pt_i : indices_groups[k]) {
			std::cout << polylines_intersections[pt_i] << ",";
			vertices_to_snapped_vertices[polylines_intersections[pt_i]] = new_v;
		}
		std::cout << "  \t moving to " << new_v << std::endl;

	}

  	return should_snap;
}

std::vector<Polyline_2> CreasePattern::merge_nearby_polylines_intersections(std::vector<Polyline_2>& polylines) {
	std::vector<Polyline_2> new_polylines;
	// First find the vertices intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(polylines.begin(), polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
  	std::cout << "initial number of intersections = " << polylines_intersections.size() << std::endl;

  	double bbox_max_len = std::max(std::abs(CGAL::to_double(bbox.xmax()-bbox.xmin())),std::abs(CGAL::to_double(bbox.ymax()-bbox.ymin())));
	Number_type dist_threshold_pow2(pow(bbox_max_len/100,2));  	
  	std::map<Point_2, Point_2> vertices_to_snapped_vertices; get_snapped_vertices_locations(polylines_intersections, dist_threshold_pow2, vertices_to_snapped_vertices);

  	Geom_traits_2::Construct_curve_2 polyline_construct = geom_traits_2.construct_curve_2_object();
  	// Add the intersection points to each curve that they are on
	// Go through segments, if the point is on the segment then split the segment, and then create 2 new segments
	// with the new snapped vertex
	int curve_i = 0;
	for (auto poly: polylines) {
		//for (auto subcurve = poly.subcurves_begin(); subcurve != poly.subcurves_end(); subcurve++) 
		std::list<Segment_2> seg_list;
		for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
			bool sing_on_segment = false;
			// This is a different type of cgal segement class, that supports "has_on" queries 
			//		(but for some reasone doesn't work with arrangements)
			CGAL::Segment_2<Kernel> tmp_seg(seg_i->source(),seg_i->target());
			//int sing_num = 0;
			//for (auto v: polylines_intersections) {
			for (int sing_num = 0; sing_num < polylines_intersections.size(); sing_num++) {
				Point_2 v = polylines_intersections[sing_num];
				// Don't split if its at the boundary of the segment
				if ((v == tmp_seg.source()) || (v == tmp_seg.target()) ) continue;
				if (tmp_seg.has_on(v)) {
					std::cout << "curve number " << curve_i << " with singularity number " << sing_num << " since vertex " << v << " is on segment " << *seg_i << std::endl; 
					sing_on_segment = true;
					
					auto snipped_v = vertices_to_snapped_vertices[v];
					seg_list.push_back(Segment_2(seg_i->source(), snipped_v));
					seg_list.push_back(Segment_2(snipped_v, seg_i->target()));
					break;
				}
			}
			if (!sing_on_segment) {
				seg_list.push_back(*seg_i);
			}
		}
		new_polylines.push_back(polyline_construct(seg_list.begin(), seg_list.end()));
		curve_i++;
	}
	return new_polylines;
}

std::vector<Polyline_2> CreasePattern::snap_nearby_polylines_start_end_starting_points(std::vector<Polyline_2>& polylines, 
							std::vector<Point_2>& intersections) {
	std::vector<std::vector<Point_2>> poly_pts;
	for (auto poly: polylines) {std::vector<Point_2> pts; PatternBoundary::polyline_to_points(poly, pts); poly_pts.push_back(pts);}
	std::vector<Number_type> x_coords(2*polylines.size()); std::vector<Number_type> y_coords(2*polylines.size());
	int cnt = 0;
	for (auto pts : poly_pts) {
		x_coords[cnt] = pts[0].x(); y_coords[cnt] = pts[0].y(); cnt++;
		x_coords[cnt] = pts.back().x(); y_coords[cnt] = pts.back().y(); cnt++;
	}
	// Save for every x point (and y point) a mapping from their number to a possibly snapped one, and then go through these points and change them
	double bbox_max_len = std::max(std::abs(CGAL::to_double(bbox.xmax()-bbox.xmin())),std::abs(CGAL::to_double(bbox.ymax()-bbox.ymin())));
	Number_type dist_threshold_pow2(pow(bbox_max_len/100,2));  	
	std::map<Number_type, Number_type> coord_to_snapped_x = snap_coords(x_coords, dist_threshold_pow2);
	std::map<Number_type, Number_type> coord_to_snapped_y = snap_coords(y_coords, dist_threshold_pow2);

	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	std::vector<Polyline_2> snapped_poly; std::vector<Point_2> endpoints;
	for (auto pts : poly_pts) {
		pts[0].x() = coord_to_snapped_x[pts[0].x()]; pts[0].y() = coord_to_snapped_y[pts[0].y()];
		pts.back().x() = coord_to_snapped_x[pts.back().x()]; pts.back().y() = coord_to_snapped_y[pts.back().y()];
		endpoints.push_back(pts[0]); endpoints.push_back(pts.back());
		snapped_poly.push_back(polyline_construct(pts.begin(), pts.end()));
	}

	// pretty inefficient but we don't have many endpoints
	for (int i = 0; i < endpoints.size(); i++) {
		for (int j = i+1; j < endpoints.size(); j++) {
			if (endpoints[i] == endpoints[j]) intersections.push_back(endpoints[i]);
		}
	}
	std::unique(intersections.begin(),intersections.end());
	return snapped_poly;
}

std::map<Number_type, Number_type> CreasePattern::snap_coords(std::vector<Number_type>& coords, Number_type threshold) {
	std::map<Number_type, Number_type> vertices_to_snapped_vertices;
	std::vector<std::set<int> > indices_groups; std::vector<bool> in_group(coords.size(), false);
  	for (int i = 0; i < coords.size(); i++) {
  		//std::cout << "vertex at: " << polylines_intersections[i] << std::endl;
  		std::set<int> nearby_points; nearby_points.insert(i);
  		for (int j = i+1; j < coords.size(); j++) {
  			auto diff = coords[i]-coords[j];
  			if ( diff*diff < threshold) {
  				nearby_points.insert(j);
  				//should_snap = true;
  			}
  		}
  		// see if j was already added somewhere
  		bool new_set = true;
  		for (int k = 0; k < indices_groups.size(); k++) {
  			if (indices_groups[k].count(i)) {
  				new_set = false;
  				for (auto idx: nearby_points) indices_groups[k].insert(idx);
  				continue;
  			}
  		}
  		if (new_set) indices_groups.push_back(nearby_points);
  	}

  	for (int k = 0; k < indices_groups.size(); k++) {
		// Compute the average of the points there
		Number_type new_v(0);
		for (auto pt_i : indices_groups[k]) {new_v += coords[pt_i];}
		new_v = new_v / Number_type(indices_groups[k].size());

		std::cout << "Vertices group " << k << ": "; 
		for (auto pt_i : indices_groups[k]) {
			std::cout << coords[pt_i] << ",";
			vertices_to_snapped_vertices[coords[pt_i]] = new_v;
		}
		std::cout << "  \t moving to " << new_v << std::endl;
	}
	return vertices_to_snapped_vertices;
}

void CreasePattern::get_visualization_mesh_and_edges(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& face_colors,
				Eigen::MatrixXd& edge_pts1, Eigen::MatrixXd& edge_pts2) {
	
	PlanarArrangement grid_with_poly(orthogonalGrid); grid_with_poly.add_polylines(initial_fold_polylines); grid_with_poly.add_polylines(initial_bnd_polylines);
	PlanarArrangement grid_with_snapped(orthogonalGrid);
	grid_with_snapped.add_polylines(clipped_fold_polylines); grid_with_snapped.add_polylines(clipped_bnd_polylines);
	
	//std::vector<PlanarArrangement*> arrangements = {&initial_arrangement, &grid_with_poly};
	//std::vector<PlanarArrangement*> arrangements = {&initial_arrangement, &grid_with_poly, &grid_with_snapped ,&clipped_grid_arrangement};
	//std::vector<PlanarArrangement*> arrangements = {&grid_with_poly, &grid_with_snapped ,&clipped_grid_arrangement};
	std::vector<PlanarArrangement*> arrangements = {&grid_with_snapped ,&clipped_grid_arrangement};
	double spacing = 1.05*CGAL::to_double(bbox.xmax()-bbox.xmin());
	get_multiple_arrangements_visualization_mesh(arrangements, spacing, V, F,face_colors);
	get_multiple_arrangements_visualization_edges(arrangements, spacing, edge_pts1, edge_pts2);
	
}

void CreasePattern::bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline) {
	Point_2 pt1(bbox.xmin(),bbox.ymin()),pt2(bbox.xmin(),bbox.ymax()),pt3(bbox.xmax(),bbox.ymax()),pt4(bbox.xmax(),bbox.ymin());
	std::list<Point_2> pts = {pt1,pt2,pt3,pt4,pt1}; // circular list
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	//std::list<Point_2> pts2 = {pt1,pt2,pt3,Point_2(0.5*(pt3.x()+pt4.x()),0.5*(pt3.x()+pt4.x())),pt4,pt1}; // circular list
	polyline = polyline_construct(pts.begin(), pts.end());
}

bool CreasePattern::is_polyline_closed_with_tolerance(const Polyline_2& poly, Number_type threshold) {
	int seg_n = poly.subcurves_end()-poly.subcurves_begin();
	auto first_pt = poly.subcurves_begin()->source(), last_pt = (poly.subcurves_begin()+(seg_n-1))->target();
	bool is_closed = (CGAL::squared_distance(first_pt,last_pt) < threshold);
	return is_closed;
}

std::vector<Point_2> CreasePattern::align_crease_vertices_x_y_with_boundary(PatternBoundary& patternBounary, 
								const std::vector<Point_2>& crease_vertices, int number_of_poly_intersections,
								 Number_type& threshold) {
	std::vector<Point_2> snapped_crease_vertices = crease_vertices; std::vector<bool> has_snapped(crease_vertices.size());
	for (int i = 0; i < snapped_crease_vertices.size(); i++) {
		auto pt = snapped_crease_vertices[i];
		std::vector<Point_2> vertical_horizontal_int = patternBounary.get_vertical_and_horizontal_intersections(pt);
		for (int j = i+1; j < snapped_crease_vertices.size(); j++) {
			if (has_snapped[j]) continue;
			bool has_snapped_j = false;
			// Get intersections from vertical and horizontal rays eminating from the pt
			if (j < number_of_poly_intersections) {
				// This is an intersection point, so snapping should be done just by setting x and y
				auto diff_x = pt.x()-snapped_crease_vertices[j].x(); auto diff_y = pt.y()-snapped_crease_vertices[j].y();
				Number_type snapped_x = snapped_crease_vertices[j].x(); Number_type snapped_y = snapped_crease_vertices[j].y();
				if (diff_x*diff_x < threshold) {snapped_x  = pt.x(); has_snapped_j = true;}
				if (diff_y*diff_y < threshold) {snapped_y  = pt.y(); has_snapped_j = true;}
				snapped_crease_vertices[j] = Point_2(snapped_x,snapped_y);
			} else {
				for (auto intersection : vertical_horizontal_int) {
					std::cout << "snapped_crease_vertices[j] = " << snapped_crease_vertices[j] << " intersection = " << intersection << " CGAL::squared_distance(pt, intersection) = " << CGAL::squared_distance(pt, intersection) << std::endl;
					if (CGAL::squared_distance(snapped_crease_vertices[j], intersection) < threshold) {
						std::cout << "snapping " << snapped_crease_vertices[j] << " to " << intersection << std::endl;
						snapped_crease_vertices[j] = intersection;
						has_snapped_j = true;
						break;
					}
				}
			}
			has_snapped[j] = has_snapped_j;
		}
		has_snapped[i] = true;
	}
	return snapped_crease_vertices;
}

std::vector<Polyline_2> CreasePattern::snap_polylines_start_end_to_vertices(std::vector<Polyline_2>& polylines, std::vector<Point_2>& vertices,
		Number_type& threshold) {
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	std::vector<Polyline_2> snapped_polylines;
	for (auto poly: polylines) {
		std::vector<Point_2> pts; PatternBoundary::polyline_to_points(poly, pts);
		for (auto v: vertices) {
			if (CGAL::squared_distance(v, pts[0]) < threshold ) pts[0] = v;
			if (CGAL::squared_distance(v, pts.back()) < threshold ) pts.back() = v;
		}
		snapped_polylines.push_back(polyline_construct(pts.begin(),pts.end()));
	}
	return snapped_polylines;
}
