#pragma once

#include <Eigen/Dense>
#include "ArrangementDefs.h"

class PlanarArrangement {
  
public:
  PlanarArrangement() : arr(&traits){};
  PlanarArrangement(const PlanarArrangement& arrangement) : arr(arrangement.arr) {}

  // multiple polylines
  void add_segments(const std::vector<Segment_2>& segments);
  void add_segment(const Segment_2& segment);
  void add_polylines(const std::vector<Polyline_2>& polylines);
  void add_polyline(const Polyline_2& polylines);

  // Used to visualize the arrangement with libigl
  void get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors);
  void get_visualization_edges(Eigen::MatrixXd& bnd_pts1, Eigen::MatrixXd& bnd_pts2);

  void get_faces_pts(std::vector<std::vector<Point_2>>& pts) const;
  void get_faces_polygons(std::vector<Polygon_2>& polygons) const;
  void get_faces_adjacency_list(std::vector<std::vector<int>>& A) const;
  void get_faces_polygons_with_holes(std::vector<Polygon_with_holes_2>& polygons) const;


  // true if it exists, false otherwise
  bool locate_point_on_vertex(const Point_2& pt, Vertex_const_handle& v);

  Arrangement_2* get_arrangement_internal() {return &arr;};

  int get_faces_n();
  int get_vertices_n();

private:
	void get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p) const;
  static void get_face_vertices_from_circulator_iter(Arrangement_2::Ccb_halfedge_const_circulator circ, 
            std::vector<Point_2>& p);
  static void get_face_vertices_from_circulator_iter2(Arrangement_2::Ccb_halfedge_const_circulator circ, 
            std::vector<Point_2>& p);

	Geom_traits_2 traits;
	Arrangement_2 arr;
};

void get_multiple_arrangements_visualization_mesh(std::vector<PlanarArrangement*> arrangements, double spacing, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
             Eigen::MatrixXd& colors);

void get_multiple_arrangements_visualization_edges(std::vector<PlanarArrangement*> arrangements, double spacing, 
  Eigen::MatrixXd& bnd_pts1, Eigen::MatrixXd& bnd_pts2);

void sort_segments(const std::vector<Segment_2>& unsorted_seg, const Point_2& firstPoint,
                std::vector<Segment_2>& sorted_seg);