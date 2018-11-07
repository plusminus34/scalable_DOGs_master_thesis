#pragma once

#include <Eigen/Dense>

#include "CreasePattern.h"
#include "../Dog/Dog.h"

Dog dog_from_crease_pattern(const CreasePattern& dogCreasePattern);
void init_grid_polygons(const CreasePattern& creasePattern,std::vector<Polygon_2>& gridPolygons);
void set_sqr_in_polygon(const CreasePattern& creasePattern, std::vector<Polygon_2>& gridPolygons, std::vector<std::vector<bool>>& sqr_in_polygon);
void generate_mesh(const CreasePattern& creasePattern, const std::vector<Polygon_2>& gridPolygons,
					 const std::vector<std::vector<bool>>& sqr_in_polygon, 
					 std::vector<Eigen::MatrixXd>& submeshVList, std::vector<Eigen::MatrixXi>& submeshFList, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void init_mesh_vertices_and_faces_from_grid(const CreasePattern& creasePattern, Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF);

void generate_constraints(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
						const std::vector<Eigen::MatrixXi>& submeshFList, DogFoldingConstraints& foldingConstraints);
Eigen::MatrixXi generate_rendered_mesh_faces(const CreasePattern& creasePattern, const Eigen::MatrixXd& V,
			 const DogFoldingConstraints& foldingConstraints);