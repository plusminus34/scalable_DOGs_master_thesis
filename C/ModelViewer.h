#pragma once

#include "ModelState.h"
//#include "Gui/DogEditor.h"
#include "DeformationController.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMeshWire = 0,
	ViewModeMesh = 1,
	ViewModeCreases = 2,
	ViewModeGauss = 3,
	CreasesSVGReader = 4,
	ViewRulings = 5,
	ViewWallpaper = 6
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState, const DeformationController& DC);

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode prevMode;
	ViewMode viewMode;
	bool render_pos_const;
	bool render_curved_folding_properties = false;

	double rulings_length = 1;
	int rulings_mod = 1;
	double rulings_planar_eps = 0.05;
	bool new_rulings = false;
	int wallpaper_res = 5;
	bool show_curves = true;
private:
	void render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer);
	void render_crease_pattern(igl::opengl::glfw::Viewer& viewer);
	void render_crease_pattern_svg_reader(igl::opengl::glfw::Viewer& viewer);
	void render_gauss_map(igl::opengl::glfw::Viewer& viewer);
	void render_positional_constraints(igl::opengl::glfw::Viewer& viewer);
	void render_edge_points_constraints(igl::opengl::glfw::Viewer& viewer);
	void render_dog_wireframe(igl::opengl::glfw::Viewer& viewer);
	void center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF);

	void render_curved_folding_normals(igl::opengl::glfw::Viewer& viewer);

	void render_mesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& Vren, const Eigen::MatrixXi& Fren);
	void render_wallpaper(igl::opengl::glfw::Viewer& viewer);
	void reset_state() {first_rendering = true;}
	void clear_edges_and_points(igl::opengl::glfw::Viewer& viewer);

	const ModelState& state;
	const DeformationController& DC;
	bool first_rendering;
	bool switched_mode;

	// Used to draw the sphere without directly using opengl routines (but just igl)
	Eigen::MatrixXd sphereV; Eigen::MatrixXi sphereF;
};
