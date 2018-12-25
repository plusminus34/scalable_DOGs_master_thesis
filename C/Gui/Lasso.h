#ifndef __ex5__Lasso__
#define __ex5__Lasso__

//#include <igl/embree/EmbreeIntersector.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/Hit.h>

#include <set>

#include "../QuadMesh/Quad.h"

class Lasso {
public:
  
  Lasso(igl::opengl::glfw::Viewer& v,
         const Eigen::MatrixXd &V_,
         const Eigen::MatrixXi &F_tri);
  ~Lasso();
  void reinit();

  int strokeAdd(int mouse_x, int mouse_y);
  int strokeAddCurve(int mouse_x, int mouse_y);
  void strokeFinish(Eigen::VectorXi &selected_vertices);
  void strokeFinishCurve(int spline_pt_number);
  int pickVertex(int mouse_x, int mouse_y);
  Edge pickEdge(int mouse_x, int mouse_y);
  void set_curve_manually(std::vector< Eigen::Matrix<double, 1,3>  >& strokePoints, std::set<int>& strokeFacesSqr, std::vector<igl::Hit>& tri_hits);
  void set_curve_manually(Eigen::MatrixXd& strokePoints, Eigen::VectorXi& faces_tri, Eigen::MatrixXd& BarycentricCoords);

  //the stroke
  std::vector< Eigen::Matrix<double, 1,3>  > strokePoints;
  std::vector< Eigen::Matrix<double, 1,3>  > splinePoints;
  std::set<int> strokeFacesSqr;
  std::vector<igl::Hit> tri_hits;
  std::vector<Edge> const_edges;
  std::vector<std::pair<double,double>> edge_coordinates;

  std::vector<int> F1; std::vector<int> F2;

  private:
  const Eigen::MatrixXd &V;
  const Eigen::MatrixXi &F;

  void get_hit_closest_edge_coordinates(igl::Hit hit, const std::vector<Edge> &edges, Edge& closest_edge, std::pair<double,double>& edge_coords);
  void get_two_face_neighbour_groups();
  void set_spline_points_from_stroke_points(int spline_pt_number);
  void updateStrokeAndHitsFromSpline();

  //igl::EmbreeIntersector ei;
  igl::opengl::glfw::Viewer &viewer;
  
  std::vector<std::vector<unsigned int> > stroke2DPoints;
  double d = -1;
};

#endif /* defined(__ex5__Lasso__) */
