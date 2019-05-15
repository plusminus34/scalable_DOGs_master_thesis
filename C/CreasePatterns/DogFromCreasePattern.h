#pragma once

#include <Eigen/Dense>

#include "CreasePattern.h"
#include "../Dog/Dog.h"

typedef std::pair<int, Polygon_2> SubmeshPoly;

Dog dog_from_crease_pattern(const CreasePattern& dogCreasePattern);
void init_grid_polygons(const CreasePattern& creasePattern,std::vector<Polygon_2>& gridPolygons);
void set_sqr_in_polygon(const CreasePattern& creasePattern, std::vector<bool>& is_polygon_hole, std::vector<Polygon_2>& gridPolygons, std::vector<std::vector<bool>>& sqr_in_polygon);
void generate_mesh(const CreasePattern& creasePattern, const std::vector<Polygon_2>& gridPolygons,
					 const std::vector<std::vector<bool>>& sqr_in_polygon, 
					 std::vector<Eigen::MatrixXd>& submeshVList, std::vector<Eigen::MatrixXi>& submeshFList, 
					 std::vector<std::vector<bool>>& submeshV_is_inner, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void init_mesh_vertices_and_faces_from_grid(const CreasePattern& creasePattern, Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF);

void generate_constraints(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
						const std::vector<Eigen::MatrixXi>& submeshFList, DogEdgeStitching& edgeStitching,
						std::vector<Point_2>& constrained_pts_non_unique,const Eigen::MatrixXd& V);

void save_submesh_bnd_edge_points(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList,
				DogEdgeStitching& edgeStitching);
//submesh_to_bnd_edge

// return a sorted list of submesh to polygons (so the polygons in the first submesh are first, then in the second, etc)
void get_faces_partitions_to_submeshes(const CreasePattern& creasePattern, std::vector<bool>& is_polygon_hole, std::vector<SubmeshPoly>& submesh_polygons);


void pt_to_edge_coordinates(const Point_2& pt, const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
				std::vector<Edge>& edge_v_indices, Number_type& t_precise, std::vector<int>& submeshes_with_pt);

bool is_closed_polyline(const Polyline_2& poly);

void polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points);
void generate_V_ren_list(Eigen::MatrixXd& V, std::vector<Eigen::MatrixXd>& submeshVList,DogEdgeStitching& eS,
	 std::vector<Eigen::MatrixXd>& V_ren_list);

void build_edge_pt_to_snapped_edge_pt_mapping(DogEdgeStitching& edgeStitching, const Eigen::MatrixXd& V);

void generate_rendered_mesh_vertices_and_faces(const CreasePattern& creasePattern, const std::vector<SubmeshPoly>& submesh_polygons,
		std::vector<Eigen::MatrixXd>& V_ren_list, DogEdgeStitching& edgeStitching,Eigen::MatrixXd& V_ren, Eigen::MatrixXi& F_ren);


bool pt_inside_polygon(const Polygon_with_holes_2& poly, const Point_2& pt);

void get_faces_polygons_without_holes(const CreasePattern& creasePattern, const std::vector<bool>& is_polygon_hole,
	 std::vector<Polygon_with_holes_2>& polygons);

std::vector<bool> submesh_is_hole(const CreasePattern& creasePattern);
Number_type bbox_diff(const CGAL::Bbox_2& bbox1, const CGAL::Bbox_2& bbox2);
Number_type bbox_max_edge(const CGAL::Bbox_2& bbox);

void polygon_mesh_to_triangle_mesh_better(const std::vector<std::vector<int> > & vF, const Eigen::MatrixXd& V_ren, Eigen::MatrixXi& F);