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
  std::vector<Point_2> sing_points; sing_points.push_back(Point_2(0.75,1.5));
  //OrthogonalGrid orthGrid(bbox,3,3); // without singularity 
  OrthogonalGrid orthGrid(bbox,3,3, sing_points);
  auto new_poly = orthGrid.single_polyline_to_segments_on_grid(polyline);
  cout << "new_poly = " << new_poly << endl;

  PlanarArrangement arrangement_with_polyline(orthGrid); 
  arrangement_with_polyline.add_polyline(polyline);

  PlanarArrangement arrangement_with_snapped_polyline(orthGrid);
  arrangement_with_snapped_polyline.add_polyline(new_poly);


  // Visualize all
  std::vector<Eigen::MatrixXd> V_list; std::vector<Eigen::MatrixXi> F_list; std::vector<Eigen::MatrixXd> F_colors_list;
  Eigen::MatrixXd V,V_grid,V_grid_with_poly,V_snapped; Eigen::MatrixXi F,F_grid,F_grid_with_poly,F_snapped; Eigen::MatrixXd grid_colors,grid_with_poly_colors,snapped_colors;

  orthGrid.get_visualization_mesh(V_grid, F_grid, grid_colors);
  arrangement_with_polyline.get_visualization_mesh(V_grid_with_poly, F_grid_with_poly, grid_with_poly_colors);
  arrangement_with_snapped_polyline.get_visualization_mesh(V_snapped, F_snapped, snapped_colors);

  double spacing = CGAL::to_double(bbox.xmax()-bbox.xmin())+1;
  V_grid_with_poly.rowwise() += Eigen::RowVector3d(1*spacing,0,0);
  V_snapped.rowwise() += Eigen::RowVector3d(2*spacing,0,0);

  V_list.push_back(V_grid); V_list.push_back(V_grid_with_poly); V_list.push_back(V_snapped);
  F_list.push_back(F_grid); F_list.push_back(F_grid_with_poly); F_list.push_back(F_snapped);
  F_colors_list.push_back(grid_colors); F_colors_list.push_back(grid_with_poly_colors);

  igl::combine(V_list,F_list, V, F);
  Eigen::MatrixXd F_colors(grid_colors.rows()+grid_with_poly_colors.rows()+snapped_colors.rows(), grid_colors.cols()); // <-- D(A.rows() + B.rows(), ...)
  F_colors << grid_colors, grid_with_poly_colors, snapped_colors;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(F_colors);
  viewer.data().show_lines = false;
  viewer.launch();
}