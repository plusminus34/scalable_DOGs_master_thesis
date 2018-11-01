#include <igl/opengl/glfw/Viewer.h>

#include <igl/combine.h>

//#include "CreasePatterns/PatternBoundary.h"
#include "CreasePatterns/DogCreasePattern.h"
#include "CreasePatterns/OrthogonalGrid.h"
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
    /*

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
  */

  //PlanarArrangement arrangement; 
  //arrangement.add_segments(segments);
/*
  std::list<Point_2> polyline_pts; //square
  polyline_pts.push_back(Point_2(0.5,0));
  polyline_pts.push_back(Point_2(1./3,0.5));
  polyline_pts.push_back(Point_2(0.5,1.1));
  polyline_pts.push_back(Point_2(0.5,1.5));
  polyline_pts.push_back(Point_2(1,1.5));
  polyline_pts.push_back(Point_2(1,1));
  polyline_pts.push_back(Point_2(1.5,0.5));
  polyline_pts.push_back(Point_2(2,0.5));
  Polyline_2 polyline = polyline_construct(polyline_pts.begin(), polyline_pts.end());

  cout << "old poly = " << polyline << endl;

  CGAL::Bbox_2 bbox(0, 0, 2, 2);
  OrthogonalGrid orthGrid(bbox,3,3);
  std::vector<Point_2> sing_points; sing_points.push_back(Point_2(0.75,1.5));
  orthGrid.add_additional_grid_points(sing_points);
  orthGrid.initialize_grid();

  auto new_poly = orthGrid.single_polyline_to_segments_on_grid(polyline);
  cout << "new_poly = " << new_poly << endl;

  PlanarArrangement arrangement_with_polyline(orthGrid); 
  arrangement_with_polyline.add_polyline(polyline);

  PlanarArrangement arrangement_with_snapped_polyline(orthGrid);
  arrangement_with_snapped_polyline.add_polyline(new_poly);
  */
  Eigen::MatrixXd V,faceColors; Eigen::MatrixXi F;
  CGAL::Bbox_2 bbox(0, 0, 2, 2);
  int x_res = 3, y_res = 3;

  std::list<Point_2> polyline_pts1; //square
  polyline_pts1.push_back(Point_2(0.5,0));
  polyline_pts1.push_back(Point_2(5./6,0.5));
  polyline_pts1.push_back(Point_2(0.5,1.1));
  polyline_pts1.push_back(Point_2(0.5,1.5));
  polyline_pts1.push_back(Point_2(1,1.5));
  polyline_pts1.push_back(Point_2(1,1));
  polyline_pts1.push_back(Point_2(1.5,0.5));
  polyline_pts1.push_back(Point_2(2,0.5));
  Polyline_2 polyline1 = polyline_construct(polyline_pts1.begin(), polyline_pts1.end());

  std::list<Point_2> polyline_pts2; //square
  polyline_pts2.push_back(Point_2(0,0));
  polyline_pts2.push_back(Point_2(2,2));
  Polyline_2 polyline2 = polyline_construct(polyline_pts2.begin(), polyline_pts2.end());

  std::vector<Polyline_2> polylines = {polyline1,polyline2};
  DogCreasePattern dogCreasePattern(bbox, polylines, x_res, y_res);
  dogCreasePattern.get_visualization_mesh(V, F, faceColors);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(faceColors);
  viewer.data().show_lines = false;
  viewer.launch();
}