#include <igl/opengl/glfw/Viewer.h>

#include "CreasePatterns/PlanarArrangement.h"

using namespace std;

int main(int argc, char *argv[])
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);
  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();
  Point_2 points1[5];
  points1[0] = Point_2(0, 0);
  points1[1] = Point_2(2, 4);
  points1[2] = Point_2(3, 0);
  points1[3] = Point_2(4, 4);
  points1[4] = Point_2(6, 0);
  Polyline_2 pi1 = polyline_construct(&points1[0], &points1[5]);
  std::list<Point_2> points2;
  points2.push_back(Point_2(1, 3));
  points2.push_back(Point_2(0, 2));
  points2.push_back(Point_2(1, 0));
  points2.push_back(Point_2(2, 1));
  points2.push_back(Point_2(3, 0));
  points2.push_back(Point_2(4, 1));
  points2.push_back(Point_2(5, 0));
  points2.push_back(Point_2(6, 2));
  points2.push_back(Point_2(5, 3));
  points2.push_back(Point_2(4, 2));
  Polyline_2 pi2 = polyline_construct(points2.begin(), points2.end());
  std::vector<Segment_2> segs;
  segs.push_back(Segment_2(Point_2(0, 2), Point_2(1, 2)));
  segs.push_back(Segment_2(Point_2(1, 2), Point_2(3, 6)));
  segs.push_back(Segment_2(Point_2(3, 6), Point_2(5, 2)));
  Polyline_2 pi3 = polyline_construct(segs.begin(), segs.end());

  PlanarArrangement arrangement;
  std::vector<Polyline_2> polylines = {pi1,pi2,pi3};
  arrangement.add_polylines(polylines);

  cout << "Number of vertices = " << arrangement.get_vertices_n() << endl;
  cout << "Number of faces = " << arrangement.get_faces_n() << endl;

  Eigen::MatrixXd V; Eigen::MatrixXi F; Eigen::VectorXd face_colors;
  arrangement.get_visualization_mesh(V, F, face_colors);
  cout << "V = " << V << endl;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
