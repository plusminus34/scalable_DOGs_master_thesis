#include "ModelViewer.h"

#include "Gui/Rendering.h"

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	//viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
	viewer.data().clear();
	if (viewMode == ViewModeMesh) {
		render_wireframe(viewer, state.dog.getV(), state.quadTop);
  		viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  		viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
	} else if (viewMode == ViewModeCreases) {
		viewer.data().set_mesh(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
		viewer.data().set_face_based(true);
  		viewer.data().set_colors(state.creasesVisualization.faceColors);
  		viewer.data().add_edges(state.creasesVisualization.edge_pts1,state.creasesVisualization.edge_pts2,Eigen::RowVector3d(0,0,0));
  		viewer.data().add_edges(state.creasesVisualization.meshE1,state.creasesVisualization.meshE2,Eigen::RowVector3d(0,0,0));
  		viewer.core.align_camera_center(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
	}
}