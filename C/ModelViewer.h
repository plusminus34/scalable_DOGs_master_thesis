#pragma once

#include "ModelState.h"
//#include "Gui/DogEditor.h"
#include "DeformationController.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMesh = 0,
	ViewModeMeshWire = 1,
	ViewModeCreases = 1,
	ViewModeGauss = 2,
	CreasesSVGReader = 3,
	ViewRulings = 4
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

	void render_mesh(igl::opengl::glfw::Viewer& viewer, const Dog& dog);

	const ModelState& state;
	const DeformationController& DC;
	bool first_rendering;
	bool switched_mode;

	// Used to draw the sphere without directly using opengl routines (but just igl)
	Eigen::MatrixXd sphereV; Eigen::MatrixXi sphereF;
};
