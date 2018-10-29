#pragma once

#include <Eigen/Dense>
#include "ArrangementDefs.h"

class PlanarArrangement {
  
public:
  PlanarArrangement() : arr(&traits){};

  // multiple polylines
  void add_polylines(std::vector<Polyline_2>& polylines);
  void add_polyline(Polyline_2& polylines);

  void get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors);

  int get_faces_n();
  int get_vertices_n();

private:
	void get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p);

	Geom_traits_2 traits;
	Arrangement_2 arr;
};