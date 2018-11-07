#pragma once

#include <Eigen/Dense>

#include "PlanarArrangement.h"
#include "OrthogonalGrid.h"


class CreasePattern {
  
public:
	CreasePattern(const CGAL::Bbox_2& bbox, std::vector<Polyline_2> polylines, int x_res, int y_res, bool snap_rounding = false);
	CreasePattern(const CreasePattern& CreasePattern);
	// TODO constructor from polygon (requiring to support removal of faces from the grid)

	const PlanarArrangement& get_clipped_arrangement() const {return clipped_grid_arrangement;}
	const OrthogonalGrid& get_orthogonal_grid() const {return orthogonalGrid;}
	//const PlanarArrangement& get_grid_with_snapped() const {return grid_with_snapped;}
	const std::vector<Polyline_2>& get_clipped_polylines() const {return clipped_polylines;}

	void get_visualization_mesh_and_edges(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors,
			Eigen::MatrixXd& edge_pts1, Eigen::MatrixXd& edge_pts2);
private:
	// snap rounding (and possibly later project initial curves to boundary)
	void init_initial_arrangement_and_polylines(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, bool snap_rounding);
	void bbox_to_polyline(const CGAL::Bbox_2& bbox, Polyline_2& polyline);


	std::vector<Polyline_2> initial_polylines; // The boundary is the first one
	PlanarArrangement initial_arrangement;
	OrthogonalGrid orthogonalGrid;
	std::vector<Polyline_2> clipped_polylines;
	//PlanarArrangement grid_with_snapped;
	PlanarArrangement clipped_grid_arrangement;
	const CGAL::Bbox_2 bbox;
};