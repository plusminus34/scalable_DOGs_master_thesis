#include <igl/opengl/glfw/Viewer.h>

//#include "CreasePatterns/PatternBoundary.h"
#include "CreasePatterns/DogCreasePattern.h"
#include "CreasePatterns/PlanarArrangement.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

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

  std::vector<Polyline_2> polylines = {pi1,pi2,pi3,pi4};
  //OrthogonalGrid orthGrid; orthGrid.polylines_to_segments_on_grid(polylines);

  std::list<Point_2> points5; points5.push_back(Point_2(0.5,0)); points5.push_back(Point_2(0.5,1));
  std::list<Point_2> points6; points6.push_back(Point_2(0,0.5)); points6.push_back(Point_2(1,0.5));
  Polyline_2 pi5 = polyline_construct(points5.begin(), points5.end());
  Polyline_2 pi6 = polyline_construct(points6.begin(), points6.end());

  std::list<Point_2> points7;
  points7.push_back(Point_2(0.5, 0));
  points7.push_back(Point_2(0.5, 0.25));
  points7.push_back(Point_2(0.52, 0.7));
  //points7.push_back(Point_2(0.25, 0.25));
  points7.push_back(Point_2(1, 1));
  Polyline_2 pi7 = polyline_construct(points7.begin(), points7.end());

  std::list<Point_2> points8;
  points8.push_back(Point_2(0, 0.5));
  points8.push_back(Point_2(1, 0.5));
  //points7.push_back(Point_2(0.25, 0.25));
  //points8.push_back(Point_2(1, 0.4));
  Polyline_2 pi8 = polyline_construct(points8.begin(), points8.end());

  //std::vector<Polyline_2> grid_and_poly = {pi1,pi5,pi6,pi7};
  //std::vector<Polyline_2> grid_and_poly = {pi1,pi5,pi7};
  std::vector<Polyline_2> grid_and_poly = {pi7,pi8};

  Geom_traits_2 geom_traits_2;
  std::vector<Polyline_2_Monotone> sub_polylines;
  CGAL::compute_subcurves(grid_and_poly.begin(), grid_and_poly.end(), std::back_inserter(sub_polylines), false, geom_traits_2);
  cout << "sub_polylines.size() = " << sub_polylines.size() << endl;
  int cnt = 1;
  //for (auto mono_poly : sub_polylines) {
  for (std::vector<Polyline_2_Monotone>::iterator scv_iter = sub_polylines.begin(); scv_iter != sub_polylines.end(); scv_iter++) {
    cout << "poly line number " << cnt;
    for (auto it = scv_iter->subcurves_begin(); it != scv_iter->subcurves_end(); it++) {
        cout << " segment " << *it << ",";
    }
    cout << endl;
    cnt++;
  }

  PlanarArrangement arrangement;
  arrangement.add_polylines(grid_and_poly);

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