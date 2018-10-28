#pragma once

#include <Eigen/Dense>

class PlanarArrangement {
  
public:
  PlanarArrangement(){};

  // multiple polylines
  void add_polylines();
  void add_polyline();

  void get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd& colors);
};