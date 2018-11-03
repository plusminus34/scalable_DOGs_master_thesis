#include "DogCreasePattern.h"

#include "OrthogonalGrid.h"

#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

#include <igl/combine.h>

DogCreasePattern::DogCreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, int x_res, int y_res, bool snap_rounding) :
											bbox(bbox), orthogonalGrid(bbox, x_res, y_res) {

	std::cout << "here" << std::endl;
	init_initial_arrangement_and_polylines(bbox, polylines, snap_rounding);
	/*
	std::cout << "there" << std::endl;
	// get poly lines intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(initial_polylines.begin(), initial_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);
  	std::cout << "ok" << std::endl;
  	// Create an orthogonal grid with singularities
  	orthogonalGrid.add_additional_grid_points(polylines_intersections);
  	std::cout << "wassup" << std::endl;
  	orthogonalGrid.initialize_grid();
  	std::cout << "ok2" << std::endl;
  	
  	
	// get new polylines
	//for (int j = 1; j < initial_polylines.size(); j++) clipped_polylines.push_back(initial_polylines[j]); // don't copy the border polygon
	for (auto poly = initial_polylines.begin()+1; poly != initial_polylines.end(); poly++) {
		std::cout << "*poly = " << *poly << std::endl;
		clipped_polylines.push_back(orthogonalGrid.single_polyline_to_segments_on_grid(*poly));
	}
	std::cout << "hmm" << std::endl;
	clipped_grid_arrangement.add_polyline(initial_polylines[0]); // add the border polygon (no need to call "clip on that")
	std::cout << "hmm2" << std::endl;
	clipped_grid_arrangement.add_polylines(clipped_polylines);
	std::cout << "hmm3" << std::endl;
	*/
}

void DogCreasePattern::init_initial_arrangement_and_polylines(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, bool snap_rounding) {
	// Create an arrangement with a boundary box polygon and the polylings
	Polyline_2 boundary_poly; bbox_to_polyline(bbox, boundary_poly);
	std::cout << "boundary_poly = " << boundary_poly << std::endl;
	// init tmp polylines with the boundary poly and tmp polylines
	std::vector<Polyline_2> tmp_polylines; tmp_polylines.push_back(boundary_poly); for (auto p: polylines) tmp_polylines.push_back(p);
	if (!snap_rounding) {
		initial_polylines = tmp_polylines;
	} else {
		// Call snap rounding with a small pixel size
		auto pixel_size = 1e-8*CGAL::min(bbox.xmax()-bbox.xmin(),bbox.ymax()-bbox.ymax());
		/*
		CGAL::snap_rounding_2<Snap_Traits,Polyline_list_2,Polyline_list_2>(tmp_polylines.begin(), tmp_polylines.end(), initial_polylines, 
																			pixel_size, false, false, 5);	
																			*/
	}
	
	// Set up the initial arrangement
	std::cout << "adding " << initial_polylines.size() << " polylines" << std::endl;
	initial_arrangement.add_polylines(initial_polylines);
}

void DogCreasePattern::get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& face_colors) {
	PlanarArrangement grid_with_snapped(orthogonalGrid);
	//PlanarArrangement grid_with_poly(orthogonalGrid); grid_with_poly.add_polylines(initial_polylines);

	//grid_with_snapped.add_polylines(clipped_polylines);
	std::vector<PlanarArrangement*> arrangements = {&initial_arrangement};
	double spacing = CGAL::to_double(bbox.xmax()-bbox.xmin())+1;
	get_multiple_arrangements_visualization_mesh(arrangements, spacing, V, F,face_colors);
}

void DogCreasePattern::bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline) {
	Point_2 pt1(bbox.xmin(),bbox.ymin()),pt2(bbox.xmin(),bbox.ymax()),pt3(bbox.xmax(),bbox.ymax()),pt4(bbox.xmax(),bbox.ymin());
	std::list<Point_2> pts = {pt1,pt2,pt3,pt4,pt1}; // circular list
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	polyline = polyline_construct(pts.begin(), pts.end());
}