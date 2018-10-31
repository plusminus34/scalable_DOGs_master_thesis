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
  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();
  /*
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);
  

  std::list<Point_2> points1; //square
  points1.push_back(Point_2(0, 0));
  points1.push_back(Point_2(1, 0));
  points1.push_back(Point_2(1, 1));
  points1.push_back(Point_2(0, 1));
  points1.push_back(Point_2(0, 0));
  Polyline_2 pi1 = polyline_construct(points1.begin(), points1.end());
  */

  // Create a 2x2 grid
  std::vector<Segment_2> segments;
  // x axis
  segments.push_back(Segment_2(Point_2(0, 0), Point_2(2, 0)));
  segments.push_back(Segment_2(Point_2(0, 1), Point_2(2, 1)));
  segments.push_back(Segment_2(Point_2(0, 2), Point_2(2, 2)));

  // y axis
  segments.push_back(Segment_2(Point_2(0, 0), Point_2(0, 2)));
  segments.push_back(Segment_2(Point_2(1, 0), Point_2(1, 2)));
  segments.push_back(Segment_2(Point_2(2, 0), Point_2(2, 2)));

  PlanarArrangement arrangement; 
  arrangement.add_segments(segments);

  std::list<Point_2> polyline_pts; //square
  polyline_pts.push_back(Point_2(0.5,0));
  polyline_pts.push_back(Point_2(1./3,0.5));
  polyline_pts.push_back(Point_2(0.5,1));
  polyline_pts.push_back(Point_2(0.5,1.5));
  polyline_pts.push_back(Point_2(1,1.5));
  polyline_pts.push_back(Point_2(1,1));
  polyline_pts.push_back(Point_2(1.5,0.5));
  polyline_pts.push_back(Point_2(2,0.5));
  Polyline_2 polyline = polyline_construct(polyline_pts.begin(), polyline_pts.end());
  arrangement.add_polyline(polyline);


  Eigen::MatrixXd V; Eigen::MatrixXi F; Eigen::MatrixXd face_colors;
  arrangement.get_visualization_mesh(V, F, face_colors);


  // query a point on the polyline
  //Point_location   pl(arr);
  //Point_2 first_point = Point_2(0.5,0);
  //Point_2          p1(4, 6);
  //point_location_query(pl, p1);

  /*
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

  arrangement.get_visualization_mesh(V, F, face_colors);
  */


  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(face_colors);
  viewer.data().show_lines = false;
  viewer.launch();
}