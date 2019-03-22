#include "Rendering.h"

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2,
    bool display_border) {
  /*
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

*/
  int e_num = quadTop.E.rows();
  int e_disp_n;
  if (display_border) {
    e_disp_n = e_num;
  } else {
    e_disp_n = 0;
    for (int i = 0; i < e_num; i++) {
      int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
      if (!quadTop.is_bnd_v[v1] && !quadTop.is_bnd_v[v2]) {
        e_disp_n++;
      }
    }
  }
  E1.resize(e_disp_n,3); E2.resize(e_disp_n,3);
  int c = 0;
  for (int i = 0; i < e_num; i++) {
    int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
    if ((display_border) || (!quadTop.is_bnd_v[v1] && !quadTop.is_bnd_v[v2]) ) {
      E1.row(c) = V.row(v1); E2.row(c) = V.row(v2);
      c++;
    }
  }
}

void render_wireframe_boundary(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
                        Eigen::RowVector3d color) {
  int e_bnd_num = 0;
  int e_num = quadTop.E.rows();
  for (int i = 0; i < e_num; i++) {
      int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
      if (quadTop.is_bnd_v[v1] && quadTop.is_bnd_v[v2]) {
        e_bnd_num++;
      }
  }
  Eigen::MatrixXd E1(e_bnd_num,3), E2(e_bnd_num,3);
  int c = 0;
  for (int i = 0; i < e_num; i++) {
    int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
    if (quadTop.is_bnd_v[v1] && quadTop.is_bnd_v[v2]) {
      E1.row(c) = V.row(v1); E2.row(c) = V.row(v2);
      c++;
    }
  }
  viewer.data().add_edges(E1, E2, color);
}

void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
    bool display_border) {
  Eigen::MatrixXd E1, E2;
  get_wireframe_edges(V, quadTop, E1, E2, display_border);
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
}

void render_dog_stitching_curves(igl::opengl::glfw::Viewer& viewer, const Dog& dog, Eigen::RowVector3d color) {
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
    viewer.data().add_edges(E1, E2, color);
  }
}

void render_dog_stitching_constraints(igl::opengl::glfw::Viewer& viewer, const Dog& dog, Eigen::RowVector3d color) {
    auto eS = dog.getEdgeStitching();
    int const_n = eS.edge_const_1.size();
    const Eigen::MatrixXd& V = dog.getV();
    Eigen::MatrixXd E1(const_n,3), E2(const_n,3);
    for (int i = 0; i < const_n; i++) {
      double t = eS.edge_coordinates[i];
      EdgePoint e1(eS.edge_const_1[i], t), e2(eS.edge_const_2[i], t);
      auto pt1 = e1.t*V.row(e1.edge.v1)+(1-e1.t)*V.row(e1.edge.v2);
      auto pt2 = e2.t*V.row(e2.edge.v1)+(1-e2.t)*V.row(e2.edge.v2);

      E1.row(i) = e1.getPositionInMesh(V);
      E2.row(i) = e2.getPositionInMesh(V);
    }
    viewer.data().add_edges(E1, E2, color);
}

void plot_vertex_based_rulings(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN, 
    const QuadTopology& quadTop, double ruling_length, int rulings_mod) {

  Eigen::MatrixXd E1(0,3); Eigen::MatrixXd E2(0,3);
  get_rulings_edges(V, VN, quadTop, ruling_length, rulings_mod,E1, E2);
  //viewer.data.add_edges(E1, E2, Eigen::RowVector3d(93./255,125./255,190./255));
  //viewer.data.add_edges(E1, E2, Eigen::RowVector3d(80./255,133./255,250./255));
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0/255,0./255,0./255));
  //viewer.data.add_edges(E1, E2, Eigen::RowVector3d(135./255,206./255,250./255));
}

void get_rulings_edges(const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN, 
                        const QuadTopology& quadTop, double ruling_length, int rulings_mod,
                        Eigen::MatrixXd& E1, Eigen::MatrixXd& E2) {
    int r_num = 0;
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
    if (!quadTop.is_bnd_v[p_xf_i] && !quadTop.is_bnd_v[p_yf_i] && !quadTop.is_bnd_v[p_xb_i] && !quadTop.is_bnd_v[p_yb_i]) {
      if (p_0_i%rulings_mod == 0) {
        r_num++;
      }
    }
  }

  // plot rulings on non flat areas
  E1.resize(r_num,3);E2.resize(r_num,3);
  int c = 0; //int flat_inner_vertices = 0;
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
    if (!quadTop.is_bnd_v[p_xf_i] && !quadTop.is_bnd_v[p_yf_i] && !quadTop.is_bnd_v[p_xb_i] && !quadTop.is_bnd_v[p_yb_i]) {
      if (p_0_i%rulings_mod == 0) {
      
      Eigen::RowVector3d r;
      r = get_ruling_direction(VN, p_0_i, p_xf_i, p_xb_i, p_yf_i, p_yb_i);
      //r = get_ruling_direction(V, p_0_i, p_xf_i, p_xb_i, p_yf_i, p_yb_i);
      double scale = ruling_length*(V.row(p_xf_i)-V.row(p_0_i)).norm();
      E1.row(c) << V(p_0_i,0)-scale*r(0), V(p_0_i,1)-scale*r(1), V(p_0_i,2)-scale*r(2);
      E2.row(c) << V(p_0_i,0)+scale*r(0), V(p_0_i,1)+scale*r(1), V(p_0_i,2)+scale*r(2);
      c++;
      }
    }
  }
  /*
  // go through flat vertices and plot edges
  c = 0; E1.resize(2*flat_inner_vertices,3); E2.resize(2*flat_inner_vertices,3);
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
    if (!quadTop.is_bnd_v[p_xf_i] && !quadTop.is_bnd_v[p_yf_i] && !quadTop.is_bnd_v[p_xb_i] && !quadTop.is_bnd_v[p_yb_i]) {
   
      Eigen::RowVector3d n = VN.row(p_0_i);
      Eigen::RowVector3d nx_f = VN.row(p_xf_i);
      Eigen::RowVector3d nx_b = VN.row(p_xb_i);
      Eigen::RowVector3d ny_f = VN.row(p_yf_i);
      Eigen::RowVector3d ny_b = VN.row(p_yb_i);

      // make sure the normals are oriented correctly
      if (n.dot(nx_f) < 0) nx_f = -nx_f;
      if (n.dot(nx_b) < 0) nx_b = -nx_b;
      if (n.dot(ny_f) < 0) ny_f = -ny_f;
      if (n.dot(ny_b) < 0) ny_b = -ny_b;

      Eigen::RowVector3d nx = nx_f-nx_b;
      Eigen::RowVector3d ny = ny_f-ny_b;
      Eigen::RowVector3d r;

      double eps = 0.2;
      // if we are flat
      if ( (nx.norm() < eps) && (ny.norm() < eps) ) {
        double scale = (V.row(p_xf_i)-V.row(p_0_i)).norm();
        r = (V.row(p_xf_i) - V.row(p_xb_i)).normalized();
        E1.row(c) << V(p_0_i,0)-scale*r(0), V(p_0_i,1)-scale*r(1), V(p_0_i,2)-scale*r(2);
        E2.row(c) << V(p_0_i,0)+scale*r(0), V(p_0_i,1)+scale*r(1), V(p_0_i,2)+scale*r(2);
        c++;

        r = (V.row(p_yf_i) - V.row(p_yb_i)).normalized();
        E1.row(c) << V(p_0_i,0)-scale*r(0), V(p_0_i,1)-scale*r(1), V(p_0_i,2)-scale*r(2);
        E2.row(c) << V(p_0_i,0)+scale*r(0), V(p_0_i,1)+scale*r(1), V(p_0_i,2)+scale*r(2);
        c++;
      } 
    }
  }
  viewer.data.add_edges(E1, E2, Eigen::RowVector3d(1,0,0));
  */
}

Eigen::RowVector3d get_ruling_direction(const Eigen::MatrixXd& VN, int p_0_i, int p_xf_i, int p_xb_i, int p_yf_i, int p_yb_i) {
   Eigen::RowVector3d n = VN.row(p_0_i);
  Eigen::RowVector3d nx_f = VN.row(p_xf_i);
  Eigen::RowVector3d nx_b = VN.row(p_xb_i);
  Eigen::RowVector3d ny_f = VN.row(p_yf_i);
  Eigen::RowVector3d ny_b = VN.row(p_yb_i);


  // Old rulings
  // make sure the normals are oriented correctly
  if (nx_f.dot(nx_b) < 0) nx_b = -nx_b;
  if (ny_f.dot(ny_b) < 0) ny_b = -ny_b;

  Eigen::RowVector3d nx = nx_f-nx_b;
  Eigen::RowVector3d ny = ny_f-ny_b;
  
  /*
  if (nx_b.dot(n) < 0) nx_b = -nx_b;
  if (nx_f.dot(n) < 0) nx_f = -nx_f;
  if (ny_f.dot(n) < 0) ny_f = -ny_f;
  if (ny_b.dot(n) < 0) ny_b = -ny_b;
  // first compute n_x, and n_y
  Eigen::Vector3d nx = (nx_f-n).normalized() + (n-nx_b).normalized();
  Eigen::Vector3d ny = (ny_f-n).normalized() + (n-ny_b).normalized();
  */
  Eigen::RowVector3d r;

  double eps = 0.05;
  // if we are flat
  if ( (nx.norm() < eps) && (ny.norm() < eps) ) {
    r << 0,0,0; // don't plot any rulings
    //flat_inner_vertices++;
  } else {
    // make sure they are oriented correctly
    if (ny.dot(nx) < 0) ny = -ny;  
    r = n.cross(nx+ny).normalized(); // average cross product and then normalize
  }
  return r;
}

/*
Eigen::RowVector3d get_ruling_direction(const Eigen::MatrixXd& V, int p_0_i, int p_xf_i, int p_xb_i, int p_yf_i, int p_yb_i) {
  Eigen::RowVector3d ex_f = V.row(p_xf_i)-V.row(p_0_i);
  Eigen::RowVector3d ex_b = V.row(p_xb_i)-V.row(p_0_i);
  Eigen::RowVector3d ey_f = V.row(p_yf_i)-V.row(p_0_i);
  Eigen::RowVector3d ey_b = V.row(p_yb_i)-V.row(p_0_i);

  double cos_x = ex_f.dot(ex_b)/(ex_f.norm()*ex_b.norm());
  double cos_y = ey_f.dot(ey_b)/(ey_f.norm()*ey_b.norm());

  cos_x = min(cos_x,1.);cos_x = max(cos_x,-1.);
  cos_y = min(cos_x,1.);cos_y = max(cos_y,-1.);

  cout << "cos_x = " << cos_x << endl;
  cout << "cos_y = " << cos_y << endl;
  double alpha = acos(cos_x);
  double beta = acos(cos_y);

  cout << "alpha = " << alpha << " beta = " << beta << endl;

  double k_x = 2*sin(alpha)/(ex_f-ex_b).norm();
  double k_y = 2*sin(beta)/(ey_f-ey_b).norm();

  cout << "k_x = " << k_x << " k_y = " << k_y << endl;

  Eigen::RowVector3d r; r.setZero();
  double eps = 1e-5;
  // Planarity test
  if ((k_x > eps) || (k_y > eps)) {
    double theta = atan2(sqrt(k_y),sqrt(k_x));
    cout << " theta = " << theta << endl;

    Eigen::VectorXd t1 = (ex_f.normalized()-ex_b.normalized()).normalized();
    Eigen::VectorXd t2 = (ey_f.normalized()-ey_b.normalized()).normalized();

    //r = sin(theta)*t1+cos(theta)*t2;
    //r = -sin(theta)*t1-cos(theta)*t2;
    r = cos(M_PI/2 - theta)*t1+sin(M_PI/2 - theta)*t2;
    //r = 10*t2    ;
  }
  return r;
}
*/