#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <boost/algorithm/string.hpp>

#include <igl/pathinfo.h>

#include "Dog/Dog.h"
#include "DeformationController.h"
#include "ModelState.h"
#include "ModelViewer.h"

using namespace std;

bool is_optimizing = true;

ModelState state;
DeformationController DC;
//DogEditor dogEditor;
ModelViewer modelViewer(state, DC);

const int DEFAULT_GRID_RES = 21;
int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh 

void clear_all_and_set_default_params(igl::opengl::glfw::Viewer& viewer) {
  DC.init_from_new_dog(state.dog);
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
  //modelViewer.viewMode = ViewModeCreases;
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
  DC.single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'S':
    DC.edit_mode = DogEditor::SELECT_POSITIONAL;
    break;
  case 'D':
    DC.edit_mode = DogEditor::TRANSLATE;
    break;
  case 'Q':
    DC.apply_new_editor_constraint();
    is_optimizing = false;
    break;
  case 'Z':
    DC.reset_new_editor_constraint();
    break;
  case 'C':
    DC.reset_constraints();

    break;
  case 'R':
    DC.single_optimization();
    break;
  case 'E':
    exit(1);
    break;
  }
  return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if (modelViewer.viewMode == ViewModeMesh) return DC.dogEditor->callback_mouse_down();
  return false;
}
bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
  if (modelViewer.viewMode == ViewModeMesh) return  DC.dogEditor->callback_mouse_move(mouse_x, mouse_y);
  return false;
}
bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if (modelViewer.viewMode == ViewModeMesh) return  DC.dogEditor->callback_mouse_up();
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  if ( ((DC.has_constraints())) && is_optimizing) run_optimization();
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
    if (argc == 3) {x_res = y_res = std::stoi(argv[2]);};
    if (argc == 4) {x_res = std::stoi(argv[2]); y_res = std::stoi(argv[3]);};
    state.init_from_svg(input_path, x_res, y_res);
    //modelViewer.viewMode = ViewModeCreases;

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
  // Set up viewer
  igl::opengl::glfw::Viewer viewer;
  DC.init_viewer(viewer);
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

    // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(210, 800), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DOG", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

      // Expose an enumeration type
      ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "Mesh\0Crease pattern\0Gauss Map\0SVG Reader\0\0");
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg(viewer);
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace(viewer);
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      if (ImGui::Button("Setup curve constraints", ImVec2(-1,0))) {DC.setup_curve_constraints();is_optimizing = false;}
      
      ImGui::Combo("Edit mode", (int *)(&DC.edit_mode), "Select\0Translate\0Vertex Pairs\0Edges Angle\0Dihedral Angle\0 None\0\0");
      ImGui::Combo("Select mode", (int *)(&DC.select_mode), "Vertex Picker\0Edge point picker\0Curve picker\0\0");
      if (ImGui::Button("Apply new constraint", ImVec2(-1,0))) {DC.apply_new_editor_constraint();}
      if (ImGui::Button("Cancel new constraint", ImVec2(-1,0))) {DC.reset_new_editor_constraint();}
      ImGui::Combo("Solver type", (int *)(&DC.p.solverType), "None\0Newton Penalty\0Newton Flow\0\0");
      ImGui::InputDouble("Bending", &DC.p.bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &DC.p.isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Soft constraints", &DC.p.soft_pos_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Dihedral weight", &DC.p.dihedral_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Fold bias weight", &DC.p.fold_bias_weight, 0, 0, "%.4f");
      
      ImGui::InputDouble("Dihedral angle", &DC.dst_dihedral_angle, 0, 0, "%.4f");
      ImGui::InputDouble("Time step size", &DC.deformation_timestep_diff, 0, 0, "%.4f");
      //if (ImGui::InputDouble("Dihedral angle", &DC.fold_dihedral_angle, 0, 0, "%.4f") ) {DC.update_fold_constraints();};
      if (ImGui::InputDouble("Timestep", &DC.deformation_timestep, 0, 0, "%.4f") ) {DC.update_time_deformations();};
      ImGui::InputDouble("Merit penalty", &DC.p.merit_p);
      ImGui::InputDouble("Infeasability epsilon", &DC.p.infeasability_epsilon);
      ImGui::InputDouble("Infeasability filter", &DC.p.infeasability_filter);
      ImGui::InputInt("Max Newton iterations", &DC.p.max_newton_iters);
      ImGui::InputInt("Penalty repetitions", &DC.p.penalty_repetitions);
      //ImGui::Checkbox("Align Procrustes", &DC.p.align_procrustes);
      ImGui::Checkbox("Render curved normals", &modelViewer.render_curved_folding_properties);
      ImGui::Checkbox("Render constraints", &modelViewer.render_pos_const);
      //ImGui::InputInt("Edited component", &dogEditor.edited_mesh);

      ImGui::InputDouble("Constraints deviation", &DC.constraints_deviation);
      ImGui::InputDouble("objective", &DC.objective);
      ImGui::Checkbox("Is optimizing?", &is_optimizing);

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