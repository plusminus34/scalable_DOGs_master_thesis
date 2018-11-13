#pragma once

#include "ModelState.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMesh = 0,
	ViewModeCreases = 1
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState);

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode prevMode;
	ViewMode viewMode;
private:
	void render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer);
	void render_crease_pattern(igl::opengl::glfw::Viewer& viewer);

	const ModelState& state;
	bool first_rendering;
	bool switched_mode;
};
