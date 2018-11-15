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

void render_dog_stitching_curves(igl::opengl::glfw::Viewer& viewer, const Dog& dog) {
  const std::vector<std::vector<EdgePoint>> &stitched_curves = dog.getEdgeStitching().stitched_curves;
  const Eigen::MatrixXd& V = dog.getV();
  for (auto curve : stitched_curves) {
    int e_num = curve.size()-1;

    Eigen::MatrixXd E1(e_num,3), E2(e_num,3);
    for (int i = 0; i < e_num; i++) {
      EdgePoint e1(curve[i]), e2(curve[i+1]);
      auto pt1 = e1.t*V.row(e1.edge.v1)+(1-e1.t)*V.row(e1.edge.v2);
      auto pt2 = e2.t*V.row(e2.edge.v1)+(1-e2.t)*V.row(e2.edge.v2);

      E1.row(i) = pt1;
      E2.row(i) = pt2;
    }
    viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
  }
}