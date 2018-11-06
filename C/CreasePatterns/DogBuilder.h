#pragma once

#include <Eigen/Dense>

#include "DogCreasePattern.h"

class DogBuilder {

public:
	DogBuilder(const DogCreasePattern& dogCreasePattern);
	
private:
	void initialize();
	void init_grid_polygons();
	void set_sqr_in_polygon();
	//const CGAL::Bbox_2 init_bbox; std::vector<Polyline_2> init_polylines;
	DogCreasePattern creasePattern;
	std::vector<Polygon_2> gridPolygons;
	// Per polygon contains a flag (whether face 'i' intersects that polygon)
	std::vector<std::vector<bool>> sqr_in_polygon;
};