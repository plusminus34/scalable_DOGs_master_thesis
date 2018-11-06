#pragma once

#include <Eigen/Dense>

#include "DogCreasePattern.h"

class DogBuilder {

public:
	DogBuilder(const std::string& svgPath, int x_res, int y_res, bool snap_rounding = false);
	DogBuilder(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, int x_res, int y_res, bool snap_rounding = false);

	void get_submeshes(std::vector<Eigen::MatrixXd>& V, std::vector<Eigen::MatrixXi>& F);
	
private:
	void initialize();
	void init_grid_polygons();
	//const CGAL::Bbox_2 init_bbox; std::vector<Polyline_2> init_polylines;
	DogCreasePattern creasePattern;
	std::vector<Polygon_2> gridPolygons;
};