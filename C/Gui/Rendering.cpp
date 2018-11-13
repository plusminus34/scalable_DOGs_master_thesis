#include "Rendering.h"

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2) {
  int e_num = quadTop.E.rows();
  int e_disp_n;
  e_disp_n = e_num;

  E1.resize(e_disp_n, 3);
  E2.resize(e_disp_n, 3);
  int c = 0;
  for (int i = 0; i < e_num; i++)
  {
    int v1 = quadTop.E(i, 0), v2 = quadTop.E(i, 1);
    E1.row(c) = V.row(v1);
    E2.row(c) = V.row(v2);
    c++;
  }
}

void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop) {
  Eigen::MatrixXd E1, E2;
  get_wireframe_edges(V, quadTop, E1, E2);
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
}