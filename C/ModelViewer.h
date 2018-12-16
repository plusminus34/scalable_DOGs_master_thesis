#pragma once

#include "ModelState.h"
#include "Dog/DogSolver.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMesh = 0,
	ViewModeCreases = 1,
	ViewModeGauss = 2
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState, const DogSolver& dogSolver);

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode prevMode;
	ViewMode viewMode;
	bool render_pos_const;
private:
	void render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer);
	void render_crease_pattern(igl::opengl::glfw::Viewer& viewer);
	void render_gauss_map(igl::opengl::glfw::Viewer& viewer);
	void render_positional_constraints(igl::opengl::glfw::Viewer& viewer);
	void render_edge_points_constraints(igl::opengl::glfw::Viewer& viewer);
	void center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF);

	const ModelState& state;
	const DogSolver& solver;
	bool first_rendering;
	bool switched_mode;

	// Used to draw the sphere without directly using opengl routines (but just igl)
	Eigen::MatrixXd sphereV; Eigen::MatrixXi sphereF;
};
