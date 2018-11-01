#pragma once

#include <Eigen/Dense>
#include "ArrangementDefs.h"

class PlanarArrangement {
  
public:
  PlanarArrangement() : arr(&traits){};
  PlanarArrangement(PlanarArrangement& arrangement) : arr(arrangement.arr) {}

  // multiple polylines
  void add_segments(const std::vector<Segment_2>& segments);
  //void add_segment(const Segment_2& segment);
  void add_polylines(const std::vector<Polyline_2>& polylines);
  void add_polyline(const Polyline_2& polylines);

  // Used to visualize the arrangement with libigl
  void get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors);


  // true if it exists, false otherwise
  bool locate_point_on_vertex(const Point_2& pt, Vertex_const_handle& v);

  Arrangement_2* get_arrangement_internal() {return &arr;};

  int get_faces_n();
  int get_vertices_n();

private:
	void get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p);

	Geom_traits_2 traits;
	Arrangement_2 arr;
};

void get_multiple_arrangements_visualization_mesh(std::vector<PlanarArrangement*> arrangements, double spacing, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
             Eigen::MatrixXd& colors);