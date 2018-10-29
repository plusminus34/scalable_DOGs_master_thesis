#include <igl/opengl/glfw/Viewer.h>

#include "CreasePatterns/PlanarArrangement.h"

using namespace std;

int main(int argc, char *argv[])
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);
  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();

  std::list<Point_2> points1; //square
  points1.push_back(Point_2(0, 0));
  points1.push_back(Point_2(1, 0));
  points1.push_back(Point_2(1, 1));
  points1.push_back(Point_2(0, 1));
  points1.push_back(Point_2(0, 0));
  Polyline_2 pi1 = polyline_construct(points1.begin(), points1.end());

  std::list<Point_2> points2;
  points2.push_back(Point_2(0,0));
  //points2.push_back(Point_2(0.5,0.5));
  points2.push_back(Point_2(1,1));
  Polyline_2 pi2 = polyline_construct(points2.begin(), points2.end());


  std::list<Point_2> points3;
  points3.push_back(Point_2(1, 0));
  points3.push_back(Point_2(0, 1));
  Polyline_2 pi3 = polyline_construct(points3.begin(), points3.end());

  std::list<Point_2> points4;
  points4.push_back(Point_2(0.5, 0));
  points4.push_back(Point_2(0.5, 0.5));
  points4.push_back(Point_2(0.75, 1));
  Polyline_2 pi4 = polyline_construct(points4.begin(), points4.end());

  PlanarArrangement arrangement;
  std::vector<Polyline_2> polylines = {pi1,pi2,pi3,pi4};
  arrangement.add_polylines(polylines);

  cout << "Number of vertices = " << arrangement.get_vertices_n() << endl;
  cout << "Number of faces = " << arrangement.get_faces_n() << endl;

  Eigen::MatrixXd V; Eigen::MatrixXi F; Eigen::MatrixXd face_colors;
  arrangement.get_visualization_mesh(V, F, face_colors);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(face_colors);
  viewer.launch();
}