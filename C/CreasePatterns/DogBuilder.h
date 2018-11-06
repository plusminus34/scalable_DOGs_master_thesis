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
	void generate_mesh();

	void init_mesh_vertices_and_faces_from_grid(Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF);

	DogCreasePattern creasePattern;
	std::vector<Polygon_2> gridPolygons;
	// Per polygon contains a flag (whether face 'i' intersects that polygon)
	std::vector<std::vector<bool>> sqr_in_polygon;

	// Generated mesh
	int submesh_n;
	std::vector<Eigen::MatrixXd> submeshVList; std::vector<Eigen::MatrixXi> submeshFList;
	Eigen::MatrixXd V; Eigen::MatrixXi F;
};