#include <igl/opengl/glfw/Viewer.h>

#include <igl/combine.h>
#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/Timer.h>

#include "CreasePatterns/CreasePattern.h"
#include "CreasePatterns/OrthogonalGrid.h"
#include "CreasePatterns/PlanarArrangement.h"
#include "CreasePatterns/SVGReader.h"
#include "CreasePatterns/DogFromCreasePattern.h"

#include "Dog/Dog.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


#include "Optimization/CompositeObjective.h"
#include "Optimization/QuadraticConstraintsSumObjective.h"
#include "Optimization/CompositeConstraints.h"
#include "Optimization/PositionalConstraints.h"
#include "Dog/Objectives/DogConstraints.h"

#include "Optimization/Solvers/LBFGS.h"

using namespace std;

int main(int argc, char *argv[])
{
  Geom_traits_2 traits;
  Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();

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
  //igl::Timer timer; double t = timer.getElapsedTime();
  CreasePattern creasePattern(bbox, polylines, x_res, y_res);
  Eigen::MatrixXd edge_pts1,edge_pts2;
  Eigen::MatrixXd V_arr,faceColors; Eigen::MatrixXi F_arr;
  creasePattern.get_visualization_mesh_and_edges(V_arr, F_arr, faceColors,edge_pts1,edge_pts2);

  Dog dog(dog_from_crease_pattern(creasePattern));
  Eigen::MatrixXd dogV; Eigen::MatrixXi dogF;
  dog.get_rendering_mesh(dogV,dogF); // render the mesh
  // move it to the right
  double spacing = 3*1.05*CGAL::to_double(bbox.xmax()-bbox.xmin());
  dogV.rowwise() += Eigen::RowVector3d(spacing,0,0);

  Eigen::MatrixXi meshE_i; igl::edges(dogF,meshE_i);
  Eigen::MatrixXd meshE1; igl::slice(dogV,meshE_i.col(0),1, meshE1);
  Eigen::MatrixXd meshE2; igl::slice(dogV,meshE_i.col(1),1, meshE2);
    
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  viewer.data().set_mesh(V_arr, F_arr);
  viewer.core.align_camera_center(V_arr,F_arr);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(faceColors);


  //viewer.data().add_edges(edge1,edge2,Eigen::RowVector3d(0,0,0));
  viewer.data().add_edges(edge_pts1,edge_pts2,Eigen::RowVector3d(0,0,0));
  viewer.data().add_edges(meshE1,meshE2,Eigen::RowVector3d(0,0,0));
  viewer.data().show_lines = false;
  viewer.launch();
}