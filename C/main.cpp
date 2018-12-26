#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <boost/algorithm/string.hpp>

#include <igl/pathinfo.h>

#include "Dog/Dog.h"
#include "ModelState.h"
#include "ModelViewer.h"

using namespace std;

bool is_optimizing = true;
ModelState state;
DeformationController deformationController;
ModelViewer modelViewer(state, deformationController);

double curve_timestep_diff = 0;
double timestep = 0;
const int DEFAULT_GRID_RES = 21;

void clear_all_and_set_default_params(igl::opengl::glfw::Viewer& viewer) {
  deformationController.init_from_new_dog(viewer, state.dog, state.quadTop);
}

void save_workspace() {
  std::string filename = igl::file_dialog_save();
  if (filename.empty())
    return;
  state.save_to_workspace(filename);
}

void load_svg(igl::opengl::glfw::Viewer& viewer) {
  std::string filename = igl::file_dialog_open();
  if (filename.empty())
    return;
  int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
  state.init_from_svg(filename, x_res, y_res);
  clear_all_and_set_default_params(viewer);
  modelViewer.viewMode = ViewModeCreases;
}

void load_workspace(igl::opengl::glfw::Viewer& viewer, const std::string& path) {
  cout << "loading workspace" << endl;
  state.load_from_workspace(path);
  clear_all_and_set_default_params(viewer);
}

void load_workspace(igl::opengl::glfw::Viewer& viewer) {
  std::string filename = igl::file_dialog_open();
  if (filename.empty())
    return;
  load_workspace(viewer, filename);
}

void run_optimization() {
  /*
  if (!is_optimizing)
    return;
  
  if (curve_timestep_diff) {
    //if (dogSolver.p.curve_timestep < 0.2) {dogSolver.p.curve_timestep += curve_timestep_diff;}
    if (deformationController.curve_timestep < 1) deformationController.curve_timestep += curve_timestep_diff;
    deformationController.update_positional_constraints();
  }
  */
  deformationController.single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'S':
    deformationController.mouse_mode = Editor::SELECT;
    break;
  case 'D':
    deformationController.mouse_mode = Editor::TRANSLATE;
    break;
  case 'F':
    deformationController.single_optimization();
    viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
    break;
  case 'C':
    deformationController.reset_constraints();
    break;
  case 'E':
    exit(1);
    break;
  }
  return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if (modelViewer.viewMode == ViewModeMesh) return deformationController.callback_mouse_down();
  return false;
}
bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
  if (modelViewer.viewMode == ViewModeMesh) return  deformationController.callback_mouse_move(mouse_x, mouse_y);
  return false;
}
bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if (modelViewer.viewMode == ViewModeMesh) return  deformationController.callback_mouse_up();
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  if (deformationController.has_constraints() && is_optimizing) run_optimization();
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
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg(viewer);
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace(viewer);
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      //ImGui::Combo("Deformation type", (int *)(&deformationController.deformationType), "Dihedral Folding\0Curve\0\0");
      ImGui::Combo("Mouse mode", (int *)(&deformationController.mouse_mode), "Select\0Translate\0None\0\0");
      ImGui::Combo("Select mode", (int *)(&deformationController.select_mode), "Vertex Picker\0Path picker\0Curve picker\0\0");
      ImGui::Combo("Solver type", (int *)(&deformationController.p.solverType), "None\0Newton Penalty\0Newton Flow\0\0");
      ImGui::InputDouble("Bending", &deformationController.p.bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &deformationController.p.isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Laplacian Similarity", &deformationController.p.laplacian_similarity_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Soft constraints", &deformationController.p.soft_pos_weight, 0, 0, "%.4f");
      //if (ImGui::InputDouble("Fold angle", &deformationController.folding_angle, 0, 0, "%.4f") ) dogSolver.update_positional_constraints();
      if (ImGui::InputDouble("Curve timestep", &deformationController.curve_timestep, 0, 0, "%.4f") ) deformationController.update_positional_constraints();
      ImGui::InputDouble("Timestep diff", &curve_timestep_diff);
      ImGui::InputDouble("Merit penalty", &deformationController.p.merit_p);
      ImGui::InputDouble("Infeasability epsilon", &deformationController.p.infeasability_epsilon);
      ImGui::InputDouble("Infeasability filter", &deformationController.p.infeasability_filter);
      ImGui::InputInt("Max Newton iterations", &deformationController.p.max_newton_iters);
      ImGui::InputInt("Penalty repetitions", &deformationController.p.penalty_repetitions);
      ImGui::Checkbox("Align Procrustes", &deformationController.p.align_procrustes);
      ImGui::Checkbox("Render constraints", &modelViewer.render_pos_const);

      ImGui::InputDouble("Constraints deviation", &deformationController.constraints_deviation);
      ImGui::InputDouble("objective", &deformationController.objective);

      ImGui::End();
  };
  clear_all_and_set_default_params(viewer);
  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());

  viewer.callback_key_down = callback_key_down;
  viewer.callback_pre_draw = callback_pre_draw; // calls at each frame
  viewer.core.is_animating = true;
  viewer.core.animation_max_fps = 30;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;

  viewer.data().show_lines = false;
  viewer.launch();
}