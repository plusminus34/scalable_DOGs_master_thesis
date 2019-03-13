#include "CreasePattern.h"

#include "OrthogonalGrid.h"

#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

#include <igl/combine.h>

CreasePattern::CreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2> polylines, int x_res, int y_res, bool snap_rounding) :
											bbox(bbox), orthogonalGrid(bbox, x_res, y_res) {
	init_initial_arrangement_and_polylines(bbox, polylines, snap_rounding);
	// get poly lines intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(initial_polylines.begin(), initial_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
  	std::cout << "polylines_intersections.size() = " << polylines_intersections.size() << std::endl; 
  	for (auto p : polylines_intersections) {
  		std::cout << "Intersection at " << p << std::endl;
  	}
  	//int wait; std::cin >> wait;
  	//exit(1);
  	//for (auto pt: polylines_intersections) {std::cout << "singular pt at " << pt << std::endl;} int wait; std::cin >> wait;
  	// Create an orthogonal grid with singularities
  	orthogonalGrid.add_additional_grid_points(polylines_intersections);
  	orthogonalGrid.regularize_grid(); // make the edges as regular as possible (not always possible due to singularities)
  	orthogonalGrid.initialize_grid();
  	
  	
	// get new polylines
	//for (int j = 1; j < initial_polylines.size(); j++) clipped_polylines.push_back(initial_polylines[j]); // don't copy the border polygon
	for (auto poly = initial_polylines.begin()+1; poly != initial_polylines.end(); poly++) {
		clipped_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(*poly));
	}
	
	clipped_grid_arrangement.add_polyline(initial_polylines[0]); // add the border polygon (no need to call "clip on that")
	clipped_grid_arrangement.add_polylines(clipped_polylines);


	//grid_with_snapped = orthogonalGrid;
	//grid_with_snapped.add_polylines(clipped_polylines);
	
}

CreasePattern::CreasePattern(const CreasePattern& cP) : initial_polylines(cP.initial_polylines), initial_arrangement(cP.initial_arrangement),
					orthogonalGrid(cP.orthogonalGrid), clipped_polylines(cP.clipped_polylines), clipped_grid_arrangement(cP.clipped_grid_arrangement), 
					bbox(cP.bbox) {
	// empty
}

void CreasePattern::init_initial_arrangement_and_polylines(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, bool snap_rounding) {
	// Create an arrangement with a boundary box polygon and the polylines
	Polyline_2 boundary_poly; bbox_to_polyline(bbox, boundary_poly);
	//std::cout << "boundary_poly = " << boundary_poly << std::endl;
	// init tmp polylines with the boundary poly and tmp polylines
	std::vector<Polyline_2> tmp_polylines; tmp_polylines.push_back(boundary_poly); for (auto p: polylines) tmp_polylines.push_back(p);
	if (!snap_rounding) {
		initial_polylines = tmp_polylines;
	} else {
/*
		typedef CGAL::Quotient<CGAL::MP_Float>           Number_type2;
		typedef CGAL::Cartesian<Number_type2>             Kernel2;
		typedef CGAL::Snap_rounding_traits_2<Kernel2>     Traits2;
		typedef Kernel2::Segment_2                        Segment_22;
		typedef Kernel2::Point_2                          Point_22;
		typedef std::list<Segment_22>                     Segment_list_22;
		typedef std::list<Point_22>                       Polyline_22;
		typedef std::list<Polyline_22>                    Polyline_list_22;
		typedef CGAL::Snap_rounding_traits_2<Kernel2>     Snap_Traits2;

		std::list<Segment_22>  seg_list2;
  		Polyline_list_22 output_list;
  		seg_list2.push_back(Segment_22(Point_22(0, 0), Point_22(10, 10)));
  		seg_list2.push_back(Segment_22(Point_22(0, 10), Point_22(10, 0)));
  		seg_list2.push_back(Segment_22(Point_22(3, 0), Point_22(3, 10)));
  		seg_list2.push_back(Segment_22(Point_22(7, 0), Point_22(7, 10)));
  // Generate an iterated snap-rounding representation, where the centers of
  // the hot pixels bear their original coordinates, using 5 kd trees:
  		std::cout << "before" << std::endl;
  CGAL::snap_rounding_2<Snap_Traits2,std::list<Segment_22>::const_iterator,Polyline_list_22>
    (seg_list2.begin(), seg_list2.end(), output_list, 1.0, true, false, 5);
    std::cout << "after" << std::endl;
    	exit(1);*/

    	//typedef CGAL::Quotient<CGAL::MP_Float>           Number_type2;
		//typedef CGAL::Cartesian<Number_type2>             Kernel2;
		//typedef CGAL::Snap_rounding_traits_2<Kernel>     Traits2;
		//typedef Kernel::Segment_2                        Segment_22;
		//typedef Kernel::Point_2                          Point_22;
		//typedef std::list<Segment_2>                     Segment_list_2;
		
		//typedef CGAL::Snap_rounding_traits_2<Kernel>     Snap_Traits2;
		//typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_22;
		//typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_22;
    	
		// Call snap rounding with a small pixel size
		auto pixel_size = 1e-8*CGAL::min(bbox.xmax()-bbox.xmin(),bbox.ymax()-bbox.ymax());

		std::list<Segment_2> seg_list;
		for (auto poly: tmp_polylines) {
			//for (auto subcurve = poly.subcurves_begin(); subcurve != poly.subcurves_end(); subcurve++) 
			for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
				seg_list.push_back(*seg_i);
				//seg_list.push_back(Segment_22(CGAL::to_double(seg_i->source()),CGAL::to_double(seg_i->target())));
				/*
				Point_22 pt1(CGAL::to_double(seg_i->source().x()) ,CGAL::to_double(seg_i->source().y()) );
				Point_22 pt2(CGAL::to_double(seg_i->target().x()),CGAL::to_double(seg_i->target().y()) );
				seg_list.push_back(Segment_22(pt1,pt2));
				*/
			}
		}

		std::vector<Point_2> polylines_intersections;
		
		Geom_traits_2 geom_traits_2;
  		CGAL::compute_intersection_points(seg_list.begin(), seg_list.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
        
  		std::cout << "Segments polylines_intersections.size() = " << polylines_intersections.size() << std::endl; 
  		for (auto p : polylines_intersections) {
  			std::cout << "Intersection at " << p << std::endl;
  		}
		//for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
  		
		/**/
		// TODO: 
		// 		1) Convert initial polylines to segment lists, where you save the number of segments each one has.
		// 		2) Concatenate the segments list into a big one.
		//		3) Call snap rounding to get the same amount of output polylines.
		//		4) Create new polylines from each one of the previous ones.
		//
		//
		// First check it on easy input. Maybe 1_poly.svg. What happens to the border? 
		//	(I guess whatever snapping happens it's ok since a rectangle will remain a rectangle, but this might require updating the bbox itself.
		//	Actually not because the boundary polygon is already here and its all we need)
		
		PointVecVec out_poly_point_list;
		CGAL::snap_rounding_2<Snap_Traits,Segment_list_2::const_iterator,PointVecVec>(seg_list.begin(), seg_list.end(), out_poly_point_list, 
																			1, false, false, 5);
		polylines_intersections.clear();
		std::cout << "out_poly_point_list.size() = " << out_poly_point_list.size() << std::endl;
		Geom_traits_2::Construct_curve_2 polyline_construct = geom_traits_2.construct_curve_2_object();
		std::vector<Segment_2> snapped_seg;
		for (int i = 0; i < out_poly_point_list.size(); i++) {
			//std::cout << "out_poly_point_list[i].size() = " << out_poly_point_list[i].size() << std::endl;
			//snapped_seg.push_back(Segment_2(out_poly_point_list[i].begin(), out_poly_point_list[i].end()));
			initial_polylines.push_back(polyline_construct(out_poly_point_list[i].begin(), out_poly_point_list[i].end()));
		}
		
		CGAL::compute_intersection_points(initial_polylines.begin(), initial_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
		std::cout << "After snapped rounding: segments polylines_intersections.size() = " << polylines_intersections.size() << std::endl; 
  		for (auto p : polylines_intersections) {
  			std::cout << "Intersection at " << p << std::endl;
  		}
  		
  		std::cout << "after snapping" << std::endl;
  		//exit(1);
	}
	
	// Set up the initial arrangement
	std::cout << "adding " << initial_polylines.size() << " polylines" << std::endl;
	initial_arrangement.add_polylines(initial_polylines);
}

void CreasePattern::get_visualization_mesh_and_edges(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& face_colors,
				Eigen::MatrixXd& edge_pts1, Eigen::MatrixXd& edge_pts2) {
	PlanarArrangement grid_with_poly(orthogonalGrid); grid_with_poly.add_polylines(initial_polylines);
	PlanarArrangement grid_with_snapped(orthogonalGrid);
	grid_with_snapped.add_polylines(clipped_polylines);
	
	//std::vector<PlanarArrangement*> arrangements = {&initial_arrangement, &grid_with_poly};
	std::vector<PlanarArrangement*> arrangements = {&initial_arrangement, &grid_with_poly, &grid_with_snapped ,&clipped_grid_arrangement};
	double spacing = 1.05*CGAL::to_double(bbox.xmax()-bbox.xmin());
	get_multiple_arrangements_visualization_mesh(arrangements, spacing, V, F,face_colors);
	get_multiple_arrangements_visualization_edges(arrangements, spacing, edge_pts1, edge_pts2);
}

void CreasePattern::bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline) {
	Point_2 pt1(bbox.xmin(),bbox.ymin()),pt2(bbox.xmin(),bbox.ymax()),pt3(bbox.xmax(),bbox.ymax()),pt4(bbox.xmax(),bbox.ymin());
	std::list<Point_2> pts = {pt1,pt2,pt3,pt4,pt1}; // circular list
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	polyline = polyline_construct(pts.begin(), pts.end());
}