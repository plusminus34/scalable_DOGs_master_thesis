#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <boost/algorithm/string.hpp>


#include <igl/pathinfo.h>

#include "Dog/Dog.h"

#include "Optimization/CompositeObjective.h"
#include "Optimization/CompositeConstraints.h"
#include "Optimization/PositionalConstraints.h"
#include "Optimization/QuadraticConstraintsSumObjective.h"
#include "Optimization/Solvers/LBFGS.h"

#include "Dog/Objectives/DogConstraints.h"
#include "Dog/Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Dog/Objectives/StitchingConstraints.h"
#include "Dog/Objectives/IsometryObjective.h"
#include "Dog/Objectives/SimplifiedBendingObjective.h"
#include "Dog/Solvers/DOGFlowAndProject.h"

#include "ModelState.h"
#include "ModelViewer.h"

using namespace std;

bool is_optimizing = false;
ModelState state;
ModelViewer modelViewer(state);
DOGFlowAndProject* solver = NULL;
FoldingAnglePositionalConstraintsBuilder* angleConstraintsBuilder = NULL;

const int DEFAULT_GRID_RES = 21;
double bending_weight = 1.;
double isometry_weight = 100.;
bool fold_mesh = true;
double folding_angle = 0;
int max_lbfgs_routines = 400;
double const_obj_penalty = 100.;
int penalty_repetitions = 1;

void clear_all_and_set_default_params() {
  if (solver){delete solver;}
  solver = new DOGFlowAndProject(state.dog, 1., 1,max_lbfgs_routines,penalty_repetitions);
  if (angleConstraintsBuilder) {delete angleConstraintsBuilder;}
  const DogEdgeStitching& eS = state.dog.getEdgeStitching();
  int c_i = eS.edge_const_1.size()/2; // TODO: This logic should be inside the constraints builder..
  angleConstraintsBuilder = new FoldingAnglePositionalConstraintsBuilder(state.dog.getV(), eS);
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

void change_fold_angle() {
  angleConstraintsBuilder->set_angle(folding_angle);
  angleConstraintsBuilder->get_positional_constraints(state.b,state.bc);
}

void load_workspace() {
  std::string filename = igl::file_dialog_open();
  if (filename.empty())
    return;
  load_workspace(filename);
}

void single_optimization() {
  cout << "running a single optimization routine" << endl;
  Eigen::VectorXd x0(state.dog.getV_vector()),x;

  // Constraints
  DogConstraints dogConst(state.quadTop);

  CompositeConstraints compConst({&dogConst});

  if (state.dog.has_creases()) {
    StitchingConstraints stitchingConstraints(state.quadTop,state.dog.getEdgeStitching());
    const DogEdgeStitching& eS = state.dog.getEdgeStitching();
    compConst.add_constraints(&stitchingConstraints);

    // Check for any positional constraints (for now these will only be folding constraints)
    /*
    if (state.b.rows()) {
      PositionalConstraints posConst(state.b,state.bc);
      compConst.add_constraints(&posConst);

      FoldingAngleConstraints
    }
    */
  }

  // Objectives
  SimplifiedBendingObjective bending(state.quadTop);
  IsometryObjective isoObj(state.quadTop,x0);
  QuadraticConstraintsSumObjective constObjBesidesPos(compConst);


  //CompositeObjective compObj({&bending, &isoObj,&constObjBesidesPos}, {bending_weight,isometry_weight,const_obj_penalty});
  CompositeObjective compObj({&bending, &isoObj}, {bending_weight,isometry_weight});
  if (state.b.rows()) {
    PositionalConstraints posConst(state.b,state.bc);
    
    /*
    QuadraticConstraintsSumObjective softPosConst(posConst);
    compObj.add_objective(&softPosConst,const_obj_penalty,true);
    */
    compConst.add_constraints(&posConst);
  }
  

  //solver->solve_single_iter(x0, compObj, compConst, x);
  //solver->solve_constrained(x0, compObj, compConst, x);
  LBFGSWithPenalty lbfsgSolver(max_lbfgs_routines, penalty_repetitions);
  lbfsgSolver.solve_constrained(x0, compObj, compConst, x);

  solver->resetSmoother();
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
  modelViewer.render(viewer);
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
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DOG", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

      // Expose an enumeration type
      ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "ViewModeMesh\0ViewModeCreases\0\0");
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg();
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace();
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      ImGui::InputDouble("Bending", &bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Const obj", &const_obj_penalty, 0, 0, "%.4f");
      ImGui::Checkbox("Folding", &fold_mesh);
      if (ImGui::InputDouble("Fold angle", &folding_angle, 0, 0, "%.4f") ) change_fold_angle();
      ImGui::InputInt("Max lbfgs iter", &max_lbfgs_routines);
      ImGui::InputInt("Penalty repetitions", &penalty_repetitions);
      ImGui::Checkbox("Render constraints", &modelViewer.render_pos_const);

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