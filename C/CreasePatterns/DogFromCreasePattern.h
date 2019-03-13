#pragma once

#include <Eigen/Dense>

#include "CreasePattern.h"
#include "../Dog/Dog.h"

typedef std::pair<int, Polygon_2> SubmeshPoly;

Dog dog_from_crease_pattern(const CreasePattern& dogCreasePattern);
void init_grid_polygons(const CreasePattern& creasePattern,std::vector<Polygon_2>& gridPolygons);
void set_sqr_in_polygon(const CreasePattern& creasePattern, std::vector<Polygon_2>& gridPolygons, std::vector<std::vector<bool>>& sqr_in_polygon);
void generate_mesh(const CreasePattern& creasePattern, const std::vector<Polygon_2>& gridPolygons,
					 const std::vector<std::vector<bool>>& sqr_in_polygon, 
					 std::vector<Eigen::MatrixXd>& submeshVList, std::vector<Eigen::MatrixXi>& submeshFList, 
					 std::vector<std::vector<bool>>& submeshV_is_inner, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void init_mesh_vertices_and_faces_from_grid(const CreasePattern& creasePattern, Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF);

void generate_constraints(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
						const std::vector<Eigen::MatrixXi>& submeshFList, DogEdgeStitching& edgeStitching,
						std::vector<Point_2>& constrained_pts_non_unique,const Eigen::MatrixXd& V);

// return a sorted list of submesh to polygons (so the polygons in the first submesh are first, then in the second, etc)
void get_faces_partitions_to_submeshes(const CreasePattern& creasePattern, std::vector<SubmeshPoly>& submesh_polygons);


// The rendered mesh V_ren = [Vlist[0];Vlist[1]..;Vlist[m];V_F] where V_F are the (non unique) folding points on the edges
// Vlist is just needed to get indices inside V_ren (to know how many vertices are in)
Eigen::MatrixXi generate_rendered_mesh_faces(const CreasePattern& creasePattern, std::vector<SubmeshPoly>& submesh_polygons,
			const std::vector<Eigen::MatrixXd>& submeshVList, const Eigen::MatrixXd& V_ren, const std::vector<Point_2>& constrained_pts_non_unique);

void pt_to_edge_coordinates(const Point_2& pt, const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
				std::vector<Edge>& edge_v_indices, Number_type& t_precise);

bool is_closed_polyline(const Polyline_2& poly);

void polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points);