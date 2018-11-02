#include <igl/opengl/glfw/Viewer.h>

#include <igl/combine.h>

//#include "CreasePatterns/PatternBoundary.h"
#include "CreasePatterns/DogCreasePattern.h"
#include "CreasePatterns/OrthogonalGrid.h"
#include "CreasePatterns/PlanarArrangement.h"
#include "CreasePatterns/SVGReader.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

using namespace std;

int main(int argc, char *argv[])
{

  Eigen::MatrixXd V,faceColors; Eigen::MatrixXi F;
  CGAL::Bbox_2 bbox;
  std::vector<Polyline_2> polylines;
  int x_res = 3, y_res = 3;

  if (argc >= 2){
    string svg_path(argv[1]);
    cout << "Reading svg file " << svg_path << endl;
    read_svg_crease_pattern(svg_path, bbox, polylines);
  } else {
    
    Geom_traits_2 traits;
    Geom_traits_2::Construct_curve_2 polyline_construct =
      traits.construct_curve_2_object();

    bbox = CGAL::Bbox_2(0, 0, 2, 2);

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
    polylines = {polyline1,polyline2};
  }
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