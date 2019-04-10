#include "CreasePattern.h"

#include "OrthogonalGrid.h"

#include <igl/combine.h>

CreasePattern::CreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2> polylines, std::vector<Polyline_2> bnd_polylines,
							int x_res, int y_res) :
											bbox(bbox), orthogonalGrid(bbox, x_res, y_res) {

	// Handle polyline intersections
	initial_fold_polylines = merge_nearby_polylines_intersections(polylines);

	// Setup initial boundary (for now just a boundary box)
	//Polyline_2 boundary_poly; bbox_to_polyline(bbox, boundary_poly);
	initial_bnd_polylines = bnd_polylines;
	/*initial_fold_polylines.push_back(boundary_poly);*/ 

	// The following is just for visualization
	initial_arrangement.add_polylines(initial_fold_polylines); initial_arrangement.add_polylines(initial_bnd_polylines);

	// get poly lines intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(initial_fold_polylines.begin(), initial_fold_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
  	std::cout << "polylines_intersections.size() = " << polylines_intersections.size() << std::endl; 
  	for (auto p : polylines_intersections) {std::cout << "Intersection at " << p << std::endl;}
  	
  	//for (auto pt: polylines_intersections) {std::cout << "singular pt at " << pt << std::endl;} int wait; std::cin >> wait;
  	// Create an orthogonal grid with singularities
  	orthogonalGrid.add_additional_grid_points(polylines_intersections);
  	// add additional interesction points at the start and end of every curve, after snapping
  	//std::vector<Polyline_2> filtered_and_clipped_to_boundary_polylines
  	orthogonalGrid.initialize_grid();
  	
  	
	// Clip boundary polylines to grid
	for (auto poly = initial_bnd_polylines.begin(); poly != initial_bnd_polylines.end(); poly++) {
		bool closed_polyline = true;
		clipped_bnd_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(*poly, closed_polyline));
	}
	// Create boundary
	patternBoundary = new PatternBoundary(clipped_bnd_polylines);

	// Clip fold polylines to grid, clip and snap them to the boudnary
	double bbox_max_len = std::max(std::abs(CGAL::to_double(bbox.xmax()-bbox.xmin())),std::abs(CGAL::to_double(bbox.ymax()-bbox.ymin())));
	Number_type dist_threshold_pow2(pow(bbox_max_len/50,2));  	
	for (auto poly = initial_fold_polylines.begin(); poly != initial_fold_polylines.end(); poly++) {
		auto filtered_and_snapped = patternBoundary->filter_and_snap(*poly,dist_threshold_pow2);
		clipped_fold_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(filtered_and_snapped));
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