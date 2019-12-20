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
bool optimization_measurements = false;

ModelState state;
//DogEditor dogEditor;
ModelViewer modelViewer(state, state.DC);

const int DEFAULT_GRID_RES = 21;
int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh

void clear_all_and_set_default_params(igl::opengl::glfw::Viewer& viewer) {
  state.DC.init_from_new_dog(state.dog, state.coarse_dog, state.conversion);
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
  state.DC.single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'S':
    state.DC.edit_mode = DogEditor::SELECT_POSITIONAL;
    break;
  case 'D':
    state.DC.edit_mode = DogEditor::TRANSLATE;
    break;
  case 'Q':
    state.DC.apply_new_editor_constraint();
    is_optimizing = false;
    break;
  case 'Z':
    state.DC.reset_new_editor_constraint();
    break;
  case 'C':
    state.DC.reset_constraints();

    break;
  case 'R':
    state.DC.single_optimization();
    break;
  case 'E':
    exit(1);
    break;
  case 'U':
    if(state.DC.display_mode != display_coarse) {
      state.DC.setDog(&state.coarse_dog);
      state.DC.display_mode = display_coarse;
    } else {
      state.DC.setDog(&state.dog);
      state.DC.display_mode = display_default;
    }
    //state.DC.change_submesh();
    break;
  }
  return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if ((modelViewer.viewMode == ViewModeMesh) || (modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return state.DC.dogEditor->callback_mouse_down();
  return false;
}
bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
  if ((modelViewer.viewMode == ViewModeMesh) || (modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return  state.DC.dogEditor->callback_mouse_move(mouse_x, mouse_y);
  return false;
}
bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if ((modelViewer.viewMode == ViewModeMesh) || (modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return  state.DC.dogEditor->callback_mouse_up();
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  if ( ((state.DC.has_constraints())) && is_optimizing) run_optimization();
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
    is_optimizing = false;
    state.load_from_workspace(input_path);
    /*
    if (argc == 3) {
      optimization_measurements = true;
      state.DC.setup_optimization_measurements(argv[2]);
    }
    */

  } else if (boost::iequals(basename, "planar")) {
    int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
    if (argc > 2) {x_res = y_res = std::stoi(argv[2]);};
    state.init_from_planar(x_res,y_res);
  } else if (boost::iequals(basename, "testcase")) {
    is_optimizing = false;
    state.load_from_workspace("testcase.work");
    int num_iterations = 1000;
    if (argc > 2) num_iterations = std::stoi(argv[2]);
    state.DC.store_data(num_iterations);
  } else {
    // Assume obj/off or other types
    state.init_from_mesh(input_path);
  }
  // Set up viewer
  igl::opengl::glfw::Viewer viewer;
  state.DC.init_viewer(viewer);
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    menu.draw_viewer_menu();
    if (ImGui::Button("Load camera", ImVec2(-1,0))) viewer.load_scene();
    if (ImGui::Button("Save camera", ImVec2(-1,0))) viewer.save_scene();
  };

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
      ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "MeshWire\0Mesh\0Crease pattern\0Gauss Map\0SVG Reader\0Rulings\0Curves\0");
      if (ImGui::Button("Load svg", ImVec2(-1,0))) load_svg(viewer);
      if (ImGui::Button("Load workspace", ImVec2(-1,0))) load_workspace(viewer);
      if (ImGui::Button("Save workspace", ImVec2(-1,0))) save_workspace();
      if (ImGui::Button("Setup curve constraints", ImVec2(-1,0))) {state.DC.setup_curve_constraints();is_optimizing = false;}

      ImGui::Combo("Edit mode", (int *)(&state.DC.edit_mode), "Select\0Translate\0Vertex Pairs\0Edges Angle\0Dihedral Angle\0 MV Dihedral Angle\0None\0\0");
      ImGui::Combo("Select mode", (int *)(&state.DC.select_mode), "Vertex Picker\0Edge point picker\0Curve picker\0\0");
      //ImGui::Combo("Wallaper type", (int *)(&state.DC.wallpaperType), "XY\0XUY\0XUYU\0XYU\0");
      if (ImGui::Button("Apply new constraint", ImVec2(-1,0))) {state.DC.apply_new_editor_constraint();}
      if (ImGui::Button("Cancel new constraint", ImVec2(-1,0))) {state.DC.reset_new_editor_constraint();}
      if (ImGui::Button("Set cylindrical boundary constraints ", ImVec2(-1,0))) {state.DC.set_cylindrical_boundary_constraints();}

      ImGui::Checkbox("Z only edit", &state.DC.z_only_editing);
      ImGui::InputDouble("Bending", &state.DC.p.bending_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Isometry", &state.DC.p.isometry_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Soft constraints", &state.DC.p.soft_pos_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Dihedral weight", &state.DC.p.dihedral_weight, 0, 0, "%.4f");
//      ImGui::InputDouble("Pairs weight", &state.DC.p.pair_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Fold bias weight", &state.DC.p.fold_bias_weight, 0, 0, "%.4f");
//      ImGui::InputDouble("MV bias weight", &state.DC.p.mv_bias_weight, 0, 0, "%.4f");
      ImGui::InputDouble("Stitching weight", &state.DC.p.stitching_weight, 0, 0, "%.4f");
      ImGui::InputDouble("ADMM rho", &state.DC.p.admm_rho, 0, 0, "%.4f");
      ImGui::InputDouble("JADMM gamma", &state.DC.p.admm_gamma, 0, 0, "%.4f");
      /*
      ImGui::InputDouble("Paired boundary smoothness bending multiply", &state.DC.paired_boundary_bending_weight_mult, 0, 0, "%.4f");
      ImGui::InputDouble("Paired boundary smoothness", &state.DC.p.paired_boundary_bending_weight, 0, 0, "%.4f");
*/
      ImGui::InputDouble("Dihedral angle (src)", &state.DC.src_dihedral_angle, 0, 0, "%.4f");
      ImGui::InputDouble("Dihedral angle (dst)", &state.DC.dst_dihedral_angle, 0, 0, "%.4f");/*
      ImGui::InputInt("Curve idx", &state.DC.deformed_curve_idx);
      ImGui::InputDouble("Curve k add", &state.DC.curve_k_translation, 0, 0, "%.4f");
      ImGui::InputDouble("Curve k mult", &state.DC.curve_k_mult, 0, 0, "%.4f");
      ImGui::InputDouble("Curve t add", &state.DC.curve_t_addition, 0, 0, "%.4f");
      ImGui::InputInt("Max constrained curve points", &state.DC.max_curve_points);
*/
      ImGui::InputDouble("Time step size", &state.DC.deformation_timestep_diff, 0, 0, "%.4f");/*
      //if (ImGui::InputDouble("Dihedral angle", &state.DC.fold_dihedral_angle, 0, 0, "%.4f") ) {state.DC.update_fold_constraints();};
      if (ImGui::InputDouble("Timestep", &state.DC.deformation_timestep, 0, 0, "%.4f") ) {state.DC.update_time_deformations();};
      ImGui::InputDouble("Merit penalty", &state.DC.p.merit_p);
      ImGui::InputDouble("Infeasability epsilon", &state.DC.p.infeasability_epsilon);
      ImGui::InputDouble("Infeasability filter", &state.DC.p.infeasability_filter);
      ImGui::InputDouble("Convergence threshold", &state.DC.p.convergence_threshold);
      ImGui::InputInt("Max Newton iterations", &state.DC.p.max_newton_iters);
      //ImGui::InputInt("Penalty repetitions", &state.DC.p.penalty_repetitions);
      */
      ImGui::Checkbox("Folding mode", &state.DC.p.folding_mode);
      /*
      ImGui::Checkbox("Flip M/V sign", &state.DC.p.flip_sign);
      ImGui::Checkbox("Render creases", &modelViewer.show_curves);
      ImGui::Checkbox("Show creases oscillating plane", &modelViewer.render_curved_folding_properties);
      ImGui::Checkbox("Render constraints", &modelViewer.render_pos_const);
      //ImGui::InputInt("Edited component", &dogEditor.edited_mesh);
      */
      ImGui::Checkbox("Render conversion", &modelViewer.show_conversion);

      ImGui::InputDouble("Constraints deviation", &state.DC.constraints_deviation);
      ImGui::InputDouble("objective", &state.DC.objective);
      ImGui::Checkbox("Is optimizing?", &is_optimizing);

      if (ImGui::Checkbox("Culled view", &modelViewer.culled_view) ) {modelViewer.switched_mode = true;}
//      ImGui::InputDouble("Rulings length", &modelViewer.rulings_length);
//      ImGui::InputInt("Rulings modulo", &modelViewer.rulings_mod);
//      ImGui::InputDouble("Rulings planar threshold", &modelViewer.rulings_planar_eps);

      ImGui::Combo("Solver mode", (int *)(&state.DC.solver_mode), "Standard\0Subsolvers\0Variable Splitting ADMM\0Jacobian ADMM\0Proximal Jacobian ADMM\0Serial\0Serial 2-patch Procrustes\0Global guess (cheat)\0Global guess (coarse)\0Coarse guess + Procrustes\0Experimental (AA)\0");
      //ImGui::Checkbox("Use subsolvers", &TODOfindagoodplace); or use some dropdown menu
      if (ImGui::Button("Add test angle constraint", ImVec2(-1,0))) state.DC.add_test_angle();
      if (ImGui::Button("Add test position constraint", ImVec2(-1,0))) state.DC.add_test_position();
      if (ImGui::Button("Reset timestep", ImVec2(-1,0))) {state.DC.deformation_timestep=0;}

      ImGui::End();
  };
  clear_all_and_set_default_params(viewer);

  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  viewer.core().align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());

  viewer.callback_key_down = callback_key_down;
  viewer.callback_pre_draw = callback_pre_draw; // calls at each frame
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.data().line_width = 2;

  viewer.data().show_lines = false;
  viewer.launch();
}
