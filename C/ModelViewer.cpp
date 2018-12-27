#include "ModelViewer.h"

#include "Gui/Rendering.h"
#include <igl/slice.h>
#include <igl/per_vertex_normals.h>

using namespace std;


ModelViewer::ModelViewer(const ModelState& modelState, const DogEditor& dogEditor) : 
									state(modelState), dogEditor(dogEditor) {
	viewMode = ViewModeMesh; 
	prevMode = viewMode;
	first_rendering = true;
	render_pos_const = true;

	igl::readOBJ("../../data/sphere.obj",sphereV,sphereF);
	center_and_scale_gauss_sphere(sphereV,sphereF);
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	viewer.data().clear();
	viewer.core.background_color = Eigen::Vector4f(1, 1, 1, 1);
	switched_mode = ((viewMode != prevMode) || (first_rendering));
	prevMode = viewMode;

	if (viewMode == ViewModeMesh) {
		render_mesh_and_wireframe(viewer);
	} else if (viewMode == CreasesSVGReader) {
		render_crease_pattern_svg_reader(viewer);
	} else if (viewMode == ViewModeGauss) {
		render_gauss_map(viewer);
	}

	first_rendering = false;
}

void ModelViewer::render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer) {
	if (switched_mode) viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
	if (state.dog.has_creases()) {
		render_dog_stitching_curves(viewer, state.dog);
	} else {
		render_wireframe(viewer, state.dog.getV(), state.dog.getQuadTopology());
	}
	if (render_pos_const) {
		render_positional_constraints(viewer);
		render_edge_points_constraints(viewer);
	}

	viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
	Eigen::Vector3d diffuse; diffuse << 135./255,206./255,250./255;
    Eigen::Vector3d ambient; /*ambient = 0.05*diffuse;*/ ambient<< 0.05,0.05,0.05;
    Eigen::Vector3d specular; specular << 0,0,0;// specular << 0.1,0.1,0.1,1.;
    //viewer.data.set_colors(diffuse);
    viewer.data().uniform_colors(ambient,diffuse,specular);

    Eigen::MatrixXd VN; igl::per_vertex_normals(state.dog.getVrendering(),state.dog.getFrendering(),VN);
  	viewer.data().set_normals(VN);
}

void ModelViewer::render_crease_pattern_svg_reader(igl::opengl::glfw::Viewer& viewer) {
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
	Eigen::VectorXi b; Eigen::VectorXd bc; dogEditor.get_positional_constraints(b,bc);
	Eigen::VectorXd constrained_pts_coords_vec; igl::slice(x,b,1, constrained_pts_coords_vec);

	int pts_num = b.size()/3;
	Eigen::MatrixXd E1(pts_num,3),E2(pts_num,3);
	for (int i = 0; i < pts_num; i++) {
		E1.row(i) << constrained_pts_coords_vec(i),constrained_pts_coords_vec(pts_num+i),constrained_pts_coords_vec(2*pts_num+i);
		E2.row(i) << bc(i),bc(pts_num+i),bc(2*pts_num+i);
	}
	
	viewer.data().add_edges(E1,E2,Eigen::RowVector3d(1.,0,0));
	dogEditor.render_positional_constraints();
}

void ModelViewer::render_edge_points_constraints(igl::opengl::glfw::Viewer& viewer) {
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
	dogEditor.get_edge_point_constraints(edgePoints, edgeCoords);
	if (!edgePoints.size()) return;
	Eigen::MatrixXd currentCoords = EdgePoint::getPositionInMesh(edgePoints, state.dog.getV());
	viewer.data().add_edges(currentCoords,edgeCoords,Eigen::RowVector3d(1.,0,0));
}

void ModelViewer::render_gauss_map(igl::opengl::glfw::Viewer& viewer) {
  viewer.data().set_mesh(sphereV, sphereF);
  Eigen::Vector3d diffuse; diffuse << 0.98,0.98,0.98;
  Eigen::Vector3d ambient; ambient << 0,0,0;//0.05*diffuse;
  Eigen::Vector3d specular; specular << 0,0,0;
  viewer.data().uniform_colors(ambient,diffuse,specular);
  //viewer.core.shininess = 0;
  
  viewer.core.align_camera_center(sphereV, sphereF);
  //viewer.core.show_lines = false;

  // TODO support curved folds by looking at normal map of each one separately
  Eigen::MatrixXd VN; igl::per_vertex_normals(state.dog.getVrendering(),state.dog.getFrendering(),VN);
  //viewer.data.set_normals(VN);
  render_wireframe(viewer,VN,state.dog.getQuadTopology(), false);
  if (switched_mode) viewer.core.align_camera_center(sphereV, sphereF);
}

void ModelViewer::center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF) {
  Eigen::RowVectorXd colmean = GV.colwise().mean();
  for (int i = 0; i < GV.rows(); i++) {
    GV.row(i) = GV.row(i)-colmean; // TODO can't this be done in 1 line?
  }
  Eigen::VectorXd area_v;
  igl::doublearea(GV,GF,area_v);
  double area = area_v.sum()/2.;
  double eps = 2e-1;
  double scale = sqrt((4-eps)*M_PI/area); // make it a little bit smaller so we could see the lines
  GV = GV * scale;
}
