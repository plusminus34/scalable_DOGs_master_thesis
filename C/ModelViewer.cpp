#include "ModelViewer.h"

#include "Gui/Rendering.h"

#include <igl/slice.h>


ModelViewer::ModelViewer(const ModelState& modelState, const DogSolver& dogSolver) : state(modelState), solver(dogSolver) {
	viewMode = ViewModeMesh; 
	prevMode = viewMode;
	first_rendering = true;
	render_pos_const = true;
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	viewer.data().clear();
	viewer.core.background_color = Eigen::Vector4f(1, 1, 1, 1);
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
	if (state.dog.has_creases()) {
		render_dog_stitching_curves(viewer, state.dog);
	} else {
		render_wireframe(viewer, state.dog.getV(), state.quadTop);
	}
	if (render_pos_const) render_positional_constraints(viewer);

	viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
	Eigen::Vector3d diffuse; diffuse << 135./255,206./255,250./255;
    Eigen::Vector3d ambient; /*ambient = 0.05*diffuse;*/ ambient<< 0.05,0.05,0.05;
    Eigen::Vector3d specular; specular << 0,0,0;// specular << 0.1,0.1,0.1,1.;
    //viewer.data.set_colors(diffuse);
    viewer.data().uniform_colors(ambient,diffuse,specular);
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

void ModelViewer::render_positional_constraints(igl::opengl::glfw::Viewer& viewer) {
	Eigen::VectorXd x(state.dog.getV_vector());
	Eigen::VectorXi b; Eigen::VectorXd bc; solver.get_positional_constraints(b,bc);
	Eigen::VectorXd constrained_pts_coords_vec; igl::slice(x,b,1, constrained_pts_coords_vec);

	int pts_num = b.size()/3;
	Eigen::MatrixXd E1(pts_num,3),E2(pts_num,3);
	for (int i = 0; i < pts_num; i++) {
		E1.row(i) << constrained_pts_coords_vec(i),constrained_pts_coords_vec(pts_num+i),constrained_pts_coords_vec(2*pts_num+i);
		E2.row(i) << bc(i),bc(pts_num+i),bc(2*pts_num+i);
	}
	viewer.data().add_edges(E1,E2,Eigen::RowVector3d(1.,0,0));
}