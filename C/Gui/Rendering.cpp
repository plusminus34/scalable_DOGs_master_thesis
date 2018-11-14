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

void render_dog_boundary(igl::opengl::glfw::Viewer& viewer, const Dog& dog) {
  /*
  const std::vector<std::vector<int>>& bndLoops = dog.getVrenBoundaryLoops();
  int l_cnt = 0;
  for (auto loop : bndLoops) {
    const int loop_size = loop.size(); int c = 0;
    std::cout << "loop " << l_cnt << " with size = " << loop_size << std::endl;
    Eigen::MatrixXd E1(loop_size,3), E2(loop_size,3);
    for (int i = 0; i < loop_size; i++) {
      E1.row(c) = dog.getVrendering().row(loop[i]);
      E2.row(c) = dog.getVrendering().row(loop[(i+1)%loop_size]);
      c++;
    }
    viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
    l_cnt++;
  }
 //const std::vector<std::vector<int>>& getVrenBoundaryLoops() const {return V_ren_bnd_loops;} 
 */
}