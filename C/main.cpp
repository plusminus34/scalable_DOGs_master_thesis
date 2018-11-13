#include <igl/opengl/glfw/Viewer.h>

#include <boost/algorithm/string.hpp>

#include <igl/combine.h>
#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/Timer.h>
#include <igl/pathinfo.h>

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
#include "Optimization/Solvers/LBFGS.h"

#include "Dog/Objectives/DogConstraints.h"
#include "Dog/Objectives/IsometryObjective.h"
#include "Dog/Objectives/SimplifiedBendingObjective.h"
#include "Dog/Solvers/DOGFlowAndProject.h"

#include "ModelState.h"

using namespace std;

bool is_optimizing = false;
ModelState state;
DOGFlowAndProject* solver = NULL;

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2)
{
  int e_num = quadTop.E.rows();
  int e_disp_n;
  e_disp_n = e_num;

  E1.resize(e_disp_n, 3);
  E2.resize(e_disp_n, 3);
  int c = 0;
  for (int i = 0; i < e_num; i++)
  {
    int v1 = quadTop.E(i, 0), v2 = quadTop.E(i, 1);
    E1.row(c) = V.row(v1);
    E2.row(c) = V.row(v2);
    c++;
  }
}

void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop)
{
  Eigen::MatrixXd E1, E2;
  get_wireframe_edges(V, quadTop, E1, E2);
  viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
}

void single_optimization() {
  cout << "running a single optimization routine" << endl;
  Eigen::VectorXd x0(state.dog.getV_vector()),x;

  // Objectives
  SimplifiedBendingObjective bending(state.quadTop);
  IsometryObjective isoObj(state.quadTop,x0);
  CompositeObjective compObj({&bending, &isoObj}, {1,5});

  // Constraints
  DogConstraints dogConst(state.quadTop);

  solver->solve_single_iter(x0, compObj, dogConst, x);
  //solver->resetSmoother();
  state.dog.update_V_vector(x);
}

void run_optimization() {
  if (!is_optimizing)
    return;
  single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'F':
    single_optimization();
    viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
    break;
  }
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  run_optimization();
  render_wireframe(viewer, state.dog.getV(), state.quadTop);
  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  return false;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: dog_editor input_path" << endl << "\t with input_path pointing to a mesh, svg, or a workspace" << endl;
    exit(1);
  }
  const std::string input_path = argv[1];
  std::string dirname,basename,extension,filename; igl::pathinfo(input_path, dirname, basename, extension, filename);
  if (boost::iequals(extension, "svg")) {
    std::cout << "Reading svg " << input_path << endl;
    exit(1);
  } else if (boost::iequals(extension, "work")) {
    std::cout << "Reading workspace " << input_path << endl;
    state.load_from_workspace(input_path);
    exit(1);
  } else {
    // Assume obj/off or other types
    state.init_from_mesh(input_path);
  }
  
  solver = new DOGFlowAndProject(state.dog, 1., 1);
  /*
  // check serialization
  igl::serialize(state,"State","bla",true);
  ModelState state2;
  igl::deserialize(state2,"State","bla");
  */

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());

  viewer.callback_key_down = callback_key_down;
  viewer.callback_pre_draw = callback_pre_draw; // calls at each frame

  viewer.core.is_animating = true;
  viewer.core.animation_max_fps = 30.;

  viewer.data().show_lines = false;
  viewer.launch();
}
/*
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
*/