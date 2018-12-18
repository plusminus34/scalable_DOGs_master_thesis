#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <boost/algorithm/string.hpp>

#include <igl/pathinfo.h>

#include "Dog/Dog.h"
#include "Dog/DogSolver.h"
#include "ModelState.h"
#include "ModelViewer.h"

using namespace std;

bool is_optimizing = false;
ModelState state;
DogSolver dogSolver;
ModelViewer modelViewer(state, dogSolver);

double curve_timestep_diff = 0;
double timestep = 0;
Curve* srcCurve = NULL; Curve* targetCurve = NULL;
Eigen::MatrixXd targetCurveCoords;

const int DEFAULT_GRID_RES = 21;

void clear_all_and_set_default_params() {
  dogSolver.init_from_new_dog(state.dog, state.quadTop);
}

void save_workspace() {
  std::string filename = igl::file_dialog_save();
  if (filename.empty())
    return;
  state.save_to_workspace(filename);
}

void load_svg() {
  std::string filename = igl::file_dialog_open();
  if (filename.empty())
    return;
  int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
  state.init_from_svg(filename, x_res, y_res);
  clear_all_and_set_default_params();
  modelViewer.viewMode = ViewModeCreases;
}

void load_workspace(const std::string& path) {
  cout << "loading workspace" << endl;
  state.load_from_workspace(path);
  clear_all_and_set_default_params();
}

void load_workspace() {
  std::string filename = igl::file_dialog_open();
  if (filename.empty())
    return;
  load_workspace(filename);
}

void run_optimization() {
  if (!is_optimizing)
    return;
  
  if (curve_timestep_diff) {
    //if (dogSolver.p.curve_timestep < 0.2) {dogSolver.p.curve_timestep += curve_timestep_diff;}
    if (dogSolver.p.curve_timestep > 1) return;
    dogSolver.p.curve_timestep += curve_timestep_diff;
    dogSolver.update_positional_constraints();
  }
  dogSolver.single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'F':
    dogSolver.single_optimization();
    viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
    break;
  }
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  run_optimization();
  modelViewer.render(viewer);
  return false;
}

void draw_curve(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd curveCoords, Eigen::RowVector3d color) {
  Eigen::MatrixXd E1(curveCoords.rows()-1,3),E2(curveCoords.rows()-1,3);
  for (int i = 0; i <curveCoords.rows()-1; i++ ) {E1.row(i) = curveCoords.row(i); E2.row(i) = curveCoords.row(i+1);}
  //viewer.data().add_points(intCurveCoords,Eigen::RowVector3d(0,0,1));
  cout << "E1.rows() = " << E1.rows() << endl;
  viewer.data().add_edges(E1,E2,color);
}

bool callback_pre_draw_curve(igl::opengl::glfw::Viewer& viewer) {
  viewer.data().clear();
  if (is_optimizing) {if (timestep < 1) timestep += 0.01;}
  // interpolate src curve and draw curve with timestep

  Eigen::Matrix3d frame; frame.setIdentity(); Eigen::RowVector3d t; t.setZero();

  Curve targetCur(targetCurveCoords);
  
  std::cout << "first curve " << endl;
  targetCur.print_geometric_represenation(); 
  std::cout << "second curve " << endl;
  targetCurve->print_geometric_represenation();

  //for (auto targetCur.length)
  //int x; cout << "check " << endl; cin >> x;

  Curve intCurve(*srcCurve,targetCur, timestep);
  std::cout << "int curve " << endl;
  intCurve.print_geometric_represenation();
  Eigen::MatrixXd intCurveCoords = intCurve.getCoords(t,frame);

  draw_curve(viewer, targetCurveCoords, Eigen::RowVector3d(1,0,0));
  draw_curve(viewer, intCurveCoords, Eigen::RowVector3d(0,0,0));
  
  // draw curve
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
    int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
    if (argc > 2) {x_res = y_res = std::stoi(argv[2]);};
    state.init_from_svg(input_path, x_res, y_res);
    modelViewer.viewMode = ViewModeCreases;

  } else if (boost::iequals(extension, "work")) {
    std::cout << "Reading workspace " << input_path << endl;
    state.load_from_workspace(input_path);

  } else if (boost::iequals(basename, "curve")) { 
    std::cout << "checking curve inteprolation!" << std::endl;
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;

    int pts_n = 20;
    std:vector<double> len,curvature1,torsion1, curvature2,torsion2;
    for (int i = 0; i < pts_n-1; i++) len.push_back(1);
    for (int i = 0; i < pts_n-2; i++) {curvature1.push_back(0); curvature2.push_back(0.2);}
    //for (int i = 0; i < pts_n-3; i++) {torsion1.push_back(0); torsion2.push_back(0.05*((pts_n-3)/2.-i));}
    //for (int i = 0; i < pts_n-2; i++) {curvature1.push_back(0); curvature2.push_back(0.05*((pts_n-2)/2.-i));}
    for (int i = 0; i < pts_n-3; i++) {torsion1.push_back(0); torsion2.push_back(curvature2[i]);}
    
    srcCurve = new Curve(len,curvature1, torsion1);
    targetCurve = new Curve(len,curvature2,torsion2);//Curve targetCurve(len,curvature2, torsion2);

    Eigen::Matrix3d frame; frame.setIdentity(); Eigen::RowVector3d t; t.setZero();
    targetCurveCoords = targetCurve->getCoords(t,frame);

    viewer.callback_key_down = callback_key_down;
    viewer.callback_pre_draw = callback_pre_draw_curve; // calls at each frame
    viewer.core.is_animating = true;
    viewer.core.animation_max_fps = 30.;
    viewer.launch();

  } else if (boost::iequals(basename, "planar")) {
    int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
    if (argc > 2) {x_res = y_res = std::stoi(argv[2]);};
    state.init_from_planar(x_res,y_res);
  } else {
    // Assume obj/off or other types
    state.init_from_mesh(input_path);
  }
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

    // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(210, 500), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DOG", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

      // Expose an enumeration type
      ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "Mesh\0Crease pattern\0Gauss Map\0\0");
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg();
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace();
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      ImGui::Combo("Deformation type", (int *)(&dogSolver.p.deformationType), "Dihedral Folding\0Curve\0\0");
      ImGui::Combo("Solver type", (int *)(&dogSolver.p.solverType), "None\0ProjectFlow\0Penalty LBFGS\0LBFGS\0Newton\0\0");
      ImGui::InputDouble("Bending", &dogSolver.p.bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &dogSolver.p.isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Laplacian Similarity", &dogSolver.p.laplacian_similarity_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Diag length", &dogSolver.p.diag_length_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Const obj", &dogSolver.p.const_obj_penalty, 0, 0, "%.4f");
      if (ImGui::InputDouble("Fold angle", &dogSolver.p.folding_angle, 0, 0, "%.4f") ) dogSolver.update_positional_constraints();
      if (ImGui::InputDouble("Curve timestep", &dogSolver.p.curve_timestep, 0, 0, "%.4f") ) dogSolver.update_positional_constraints();
      ImGui::InputDouble("Timestep diff", &curve_timestep_diff);
      ImGui::InputInt("Max lbfgs iter", &dogSolver.p.max_lbfgs_routines);
      ImGui::InputInt("Penalty repetitions", &dogSolver.p.penalty_repetitions);
      ImGui::Checkbox("Project after", &dogSolver.p.project_after_flow);
      ImGui::Checkbox("Align Procrustes", &dogSolver.p.align_procrustes);
      ImGui::Checkbox("ARAP Guess", &dogSolver.p.arap_guess);
      ImGui::Checkbox("Render constraints", &modelViewer.render_pos_const);

      ImGui::InputDouble("Constraints deviation", &dogSolver.constraints_deviation);
      ImGui::InputDouble("objective", &dogSolver.objective);

      ImGui::End();
  };
  clear_all_and_set_default_params();
  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());

  viewer.callback_key_down = callback_key_down;
  viewer.callback_pre_draw = callback_pre_draw; // calls at each frame

  viewer.core.is_animating = true;
  viewer.core.animation_max_fps = 30.;

  viewer.data().show_lines = false;
  viewer.launch();
}