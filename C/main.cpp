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
  Geom_traits_2 traits;
  Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();

  Eigen::MatrixXd V,faceColors; Eigen::MatrixXi F;
  CGAL::Bbox_2 bbox;
  std::vector<Polyline_2> polylines;
  int x_res = 3, y_res = 3;

  if (argc >= 2){
    string svg_path(argv[1]);
    cout << "Reading svg file " << svg_path << endl;
    read_svg_crease_pattern(svg_path, bbox, polylines);

    if (argc > 2) {x_res = y_res = std::stoi(argv[2]);};
  } else {

    bbox = CGAL::Bbox_2(0, 0, 2, 2);

    std::list<Point_2> polyline_pts1;
    polyline_pts1.push_back(Point_2(0.5,0));
    polyline_pts1.push_back(Point_2(5./6,0.5));
    polyline_pts1.push_back(Point_2(0.5,1.1));
    polyline_pts1.push_back(Point_2(0.5,1.5));
    polyline_pts1.push_back(Point_2(1,1.5));
    polyline_pts1.push_back(Point_2(1,1));
    polyline_pts1.push_back(Point_2(1.5,0.5));
    polyline_pts1.push_back(Point_2(2,0.5));
    Polyline_2 polyline1 = polyline_construct(polyline_pts1.begin(), polyline_pts1.end());

    std::list<Point_2> polyline_pts2;
    polyline_pts2.push_back(Point_2(0,0));
    polyline_pts2.push_back(Point_2(2,2));
    Polyline_2 polyline2 = polyline_construct(polyline_pts2.begin(), polyline_pts2.end());
    polylines = {polyline1,polyline2};
  }
  DogCreasePattern dogCreasePattern(bbox, polylines, x_res, y_res);
  dogCreasePattern.get_visualization_mesh(V, F, faceColors);
  /*
  std::vector<std::vector<Point_2>> all_faces_pts;
  dogCreasePattern.initial_arrangement.get_faces_pts(all_faces_pts);
  std::vector<Point_2> bnd_pts = all_faces_pts[1];
  Eigen::MatrixXd bnd_pts1(bnd_pts.size(),3),bnd_pts2(bnd_pts.size(),3);
  
  for (int i = 0; i < bnd_pts.size();i++) {
  
    double x = CGAL::to_double(bnd_pts[i].x()), y = CGAL::to_double(bnd_pts[i].y());
    int next_i = (i+1)%bnd_pts.size();
    double next_x = CGAL::to_double(bnd_pts[next_i].x()), next_y = CGAL::to_double(bnd_pts[next_i].y());
    bnd_pts1.row(i) << x,y,0;
    bnd_pts2.row(i) << next_x,next_y,0;
  }
  */
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.core.align_camera_center(V,F);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(faceColors);
  //viewer.data().add_edges(bnd_pts1,bnd_pts2,Eigen::RowVector3d(1,0,0));
  viewer.data().show_lines = false;
  viewer.launch();
}