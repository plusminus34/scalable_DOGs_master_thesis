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
    ImGui::SetNextWindowSize(ImVec2(210, 500), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DOG", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

      // Expose an enumeration type
      ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "ViewModeMesh\0ViewModeCreases\0\0");
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg();
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace();
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      ImGui::Combo("Deformation type", (int *)(&dogSolver.p.deformationType), "Dihedral Folding\0Curve\0\0");
      ImGui::InputDouble("Bending", &dogSolver.p.bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &dogSolver.p.isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Const obj", &dogSolver.p.const_obj_penalty, 0, 0, "%.4f");
      if (ImGui::InputDouble("Fold angle", &dogSolver.p.folding_angle, 0, 0, "%.4f") ) dogSolver.update_positional_constraints();
      if (ImGui::InputDouble("Curve timestep", &dogSolver.p.curve_timestep, 0, 0, "%.4f") ) dogSolver.update_positional_constraints();
      ImGui::InputInt("Max lbfgs iter", &dogSolver.p.max_lbfgs_routines);
      ImGui::InputInt("Penalty repetitions", &dogSolver.p.penalty_repetitions);
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