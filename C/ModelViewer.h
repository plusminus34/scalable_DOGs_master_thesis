#pragma once

#include "ModelState.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMesh = 0,
	ViewModeCreases = 1
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState) : state(modelState) {viewMode = ViewModeMesh;}

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode viewMode;
private:
	const ModelState& state;
};

