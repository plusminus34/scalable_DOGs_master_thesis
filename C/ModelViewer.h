#pragma once

#include "ModelState.h"
//#include "Gui/DogEditor.h"
#include "DeformationController.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMesh = 0,
	ViewModeCreases = 1,
	ViewModeGauss = 2,
	CreasesSVGReader = 3,
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState, const DeformationController& DC);

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode prevMode;
	ViewMode viewMode;
	bool render_pos_const;
private:
	void render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer);
	void render_crease_pattern_svg_reader(igl::opengl::glfw::Viewer& viewer);
	void render_gauss_map(igl::opengl::glfw::Viewer& viewer);
	void render_positional_constraints(igl::opengl::glfw::Viewer& viewer);
	void render_edge_points_constraints(igl::opengl::glfw::Viewer& viewer);
	void center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF);

	const ModelState& state;
	const DeformationController& DC;
	bool first_rendering;
	bool switched_mode;

	// Used to draw the sphere without directly using opengl routines (but just igl)
	Eigen::MatrixXd sphereV; Eigen::MatrixXi sphereF;
};
