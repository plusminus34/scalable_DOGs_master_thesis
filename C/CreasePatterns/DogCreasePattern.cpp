#include "DogCreasePattern.h"

#include "OrthogonalGrid.h"

#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

#include <igl/combine.h>

DogCreasePattern::DogCreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, int x_res, int y_res, bool snap_rounding) :
											bbox(bbox), orthogonalGrid(bbox, x_res, y_res) {
	init_initial_arrangement_and_polylines(bbox, polylines, snap_rounding);
	// get poly lines intersections
	std::vector<Point_2> polylines_intersections;
	Geom_traits_2 geom_traits_2;
  	CGAL::compute_intersection_points(initial_polylines.begin(), initial_polylines.end(),
                                    std::back_inserter(polylines_intersections), false, geom_traits_2);

  	// Create an orthogonal grid with singularities
  	orthogonalGrid.add_additional_grid_points(polylines_intersections);
  	orthogonalGrid.initialize_grid();
  	
	// get new polylines
	//for (int j = 1; j < initial_polylines.size(); j++) clipped_polylines.push_back(initial_polylines[j]); // don't copy the border polygon
	clipped_grid_arrangement.add_polyline(initial_polylines[0]); // add the border polygon (no need to call "clip on that")
	for (auto poly = initial_polylines.begin()+1; poly != initial_polylines.end(); poly++) {
		auto other(orthogonalGrid);
		auto pts = other.single_polyline_to_segments_on_grid(*poly);

		Geom_traits_2 traits;
		Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
		Polyline_2 clipped_poly = polyline_construct(pts.begin(), pts.end());
		clipped_polylines.push_back(clipped_poly);
		//clipped_grid_arrangement.add_polyline(clipped_poly);
	}
	
	
	std::cout << "before = " << std::endl;
	//clipped_grid_arrangement.add_polylines(clipped_polylines);
	std::cout << "alive" << std::endl;
	///Arrangement_2 arr; 
	//insert(arr, clipped_polylines[0]);
	//insert(arr, clipped_polylines[1]);
	//insert(arr, clipped_polylines[2]);
	//insert(arr, initial_polylines[2]);
	//insert(arr, clipped_polylines.begin(), clipped_polylines.end());
	std::cout << "alive2 wtf" << std::endl;
	//exit(1);
}

void DogCreasePattern::init_initial_arrangement_and_polylines(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, bool snap_rounding) {
	// Create an arrangement with a boundary box polygon and the polylings
	Polyline_2 boundary_poly; bbox_to_polyline(bbox, boundary_poly);
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
	initial_arrangement.add_polylines(initial_polylines);
}

void DogCreasePattern::get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& face_colors) {
	PlanarArrangement grid_with_snapped(orthogonalGrid);
	PlanarArrangement grid_with_poly(orthogonalGrid); grid_with_poly.add_polylines(initial_polylines);

	grid_with_snapped.add_polylines(clipped_polylines);
	std::vector<PlanarArrangement*> arrangements = {&initial_arrangement,&grid_with_poly, &grid_with_snapped, &clipped_grid_arrangement};
	double spacing = CGAL::to_double(bbox.xmax()-bbox.xmin())+1;
	get_multiple_arrangements_visualization_mesh(arrangements, spacing, V, F,face_colors);

/*

	// Visualize all
	std::vector<Eigen::MatrixXd> V_list; std::vector<Eigen::MatrixXi> F_list; std::vector<Eigen::MatrixXd> F_colors_list;
	Eigen::MatrixXd V_init,V_grid,V_grid_with_poly,V_snapped; Eigen::MatrixXi F_init,F_grid,F_grid_with_poly,F_snapped; 
	Eigen::MatrixXd init_colors, grid_colors,grid_with_poly_colors,snapped_colors;

	initial_arrangement.get_visualization_mesh(V_init, F_init, init_colors);
	
/*
	arrangement_with_polyline.get_visualization_mesh(V_grid_with_poly, F_grid_with_poly, grid_with_poly_colors);
	arrangement_with_snapped_polyline.get_visualization_mesh(V_snapped, F_snapped, snapped_colors);
*/
	/*
	V_grid.rowwise() += Eigen::RowVector3d(1*spacing,0,0);
	//V_snapped.rowwise() += Eigen::RowVector3d(2*spacing,0,0);

	V_list.push_back(V_init); V_list.push_back(V_grid); //V_list.push_back(V_snapped);
	F_list.push_back(F_init); F_list.push_back(F_grid); //F_list.push_back(F_snapped);
	F_colors_list.push_back(init_colors); F_colors_list.push_back(grid_colors);

	igl::combine(V_list,F_list, V, F);
	face_colors.resize(init_colors.rows() + grid_colors.rows(), init_colors.cols()); // <-- D(A.rows() + B.rows(), ...)
	face_colors << init_colors,grid_colors;//, snapped_colors;*/
}

void DogCreasePattern::bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline) {
	Point_2 pt1(bbox.xmin(),bbox.ymin()),pt2(bbox.xmin(),bbox.ymax()),pt3(bbox.xmax(),bbox.ymax()),pt4(bbox.xmax(),bbox.ymin());
	std::list<Point_2> pts = {pt1,pt2,pt3,pt4,pt1}; // circular list
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
	polyline = polyline_construct(pts.begin(), pts.end());
}