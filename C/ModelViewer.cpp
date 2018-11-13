#include "ModelViewer.h"

#include "Gui/Rendering.h"


ModelViewer::ModelViewer(const ModelState& modelState) : state(modelState) {
	viewMode = ViewModeMesh; 
	prevMode = viewMode;
	first_rendering = true;
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	//viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
	viewer.data().clear();
	switched_mode = ((viewMode != prevMode) || (first_rendering));
	prevMode = viewMode;

	if (viewMode == ViewModeMesh) {
		render_mesh_and_wireframe(viewer);
	} else if (viewMode == ViewModeCreases) {
		render_crease_pattern(viewer);
	}

	first_rendering = false;
}

void ModelViewer::render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer) {
	render_wireframe(viewer, state.dog.getV(), state.quadTop);
	viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
	viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
}

void ModelViewer::render_crease_pattern(igl::opengl::glfw::Viewer& viewer) {
	if (state.dog.has_creases()) {
		viewer.data().set_mesh(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
		viewer.data().set_face_based(true);
		viewer.data().set_colors(state.creasesVisualization.faceColors);
		viewer.data().add_edges(state.creasesVisualization.edge_pts1,state.creasesVisualization.edge_pts2,Eigen::RowVector3d(0,0,0));
		viewer.data().add_edges(state.creasesVisualization.meshE1,state.creasesVisualization.meshE2,Eigen::RowVector3d(0,0,0));
		viewer.core.align_camera_center(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
	} else {
		if (switched_mode) std::cout << "This DOG has no creases" << std::endl;
	  	render_mesh_and_wireframe(viewer);
	}
}