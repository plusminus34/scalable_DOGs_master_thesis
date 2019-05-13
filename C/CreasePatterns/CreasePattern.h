#pragma once

#include <Eigen/Dense>

#include "PlanarArrangement.h"
#include "OrthogonalGrid.h"
#include "PatternBoundary.h"


class CreasePattern {
  
public:
	CreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2> polylines, std::vector<Polyline_2> bnd_polylines,
							 int x_res, int y_res);
	CreasePattern(const CreasePattern& CreasePattern);
	// TODO constructor from polygon (requiring to support removal of faces from the grid)

	const PlanarArrangement& get_clipped_arrangement() const {return clipped_grid_arrangement;}
	const OrthogonalGrid& get_orthogonal_grid() const {return orthogonalGrid;}
	//const PlanarArrangement& get_grid_with_snapped() const {return grid_with_snapped;}
	const std::vector<Polyline_2>& get_clipped_fold_polylines() const {return clipped_fold_polylines;}
	const std::vector<Polyline_2>& get_clipped_bnd_polylines() const {return clipped_bnd_polylines;}

	void get_visualization_mesh_and_edges(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors,
			Eigen::MatrixXd& edge_pts1, Eigen::MatrixXd& edge_pts2);

	//void get_submeshes_faces_polygons(std::vector<Polygon_2>& polygons) const;
	void get_submeshes_faces_polygons(std::vector<Polygon_with_holes_2>& polygons) const {
		get_clipped_arrangement().get_faces_polygons_with_holes(polygons);
	}

	const PatternBoundary* boundary() const {return patternBoundary;}
private:
	// snap rounding (and possibly later project initial curves to boundary)
	bool get_snapped_vertices_locations(const std::vector<Point_2>& polylines_int, Number_type threshold, std::map<Point_2, Point_2>& vertices_to_snapped_vertices);
	void bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline);
	bool is_polyline_closed_with_tolerance(const Polyline_2& polyline, Number_type threshold);

	std::vector<Polyline_2> merge_nearby_polylines_intersections(std::vector<Polyline_2>& polylines);
	std::vector<Polyline_2> snap_nearby_polylines_start_end_starting_points(std::vector<Polyline_2>& polylines, std::vector<Point_2>& intersections);
	std::vector<Polyline_2> snap_polylines_start_end_to_vertices(std::vector<Polyline_2>& polylines, std::vector<Point_2>& vertices, Number_type& threshold);
	std::vector<Polyline_2> snap_polylines_to_vertices(std::vector<Polyline_2>& polylines, std::vector<Point_2>& vertices,
		Number_type& threshold);
	std::map<Number_type, Number_type> snap_coords(std::vector<Number_type>& coords, Number_type threshold);
	std::vector<Polyline_2> snap_and_split_curves_to_starting_points(std::vector<Polyline_2>& polylines, Number_type threshold);

	std::vector<Point_2> align_crease_vertices_x_y_with_boundary(PatternBoundary& patternBounary, 
								const std::vector<Point_2>& crease_vertices, int number_of_poly_int, Number_type& threshold);


	std::vector<Polyline_2> initial_fold_polylines;
	std::vector<Polyline_2> initial_bnd_polylines;
	PlanarArrangement initial_arrangement;
	OrthogonalGrid orthogonalGrid;
	std::vector<Polyline_2> clipped_fold_polylines;
	std::vector<Polyline_2> clipped_bnd_polylines;
	//PlanarArrangement grid_with_snapped;
	PlanarArrangement clipped_grid_arrangement;
	const CGAL::Bbox_2 bbox;
	PatternBoundary* patternBoundary;
	std::vector<Point_2> submeshes_mass_center;
};