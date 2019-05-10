#include "ModelViewer.h"

#include "Gui/Rendering.h"
#include <igl/combine.h>
#include <igl/slice.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>

using namespace std;


ModelViewer::ModelViewer(const ModelState& modelState, const DeformationController& DC) : 
									state(modelState), DC(DC) {
	viewMode = ViewModeMeshWire; 
	prevMode = viewMode;
	first_rendering = true;
	render_pos_const = true;

	igl::readOBJ("../../data/sphere.obj",sphereV,sphereF);
	center_and_scale_gauss_sphere(sphereV,sphereF);
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = DC.getEditedSubmesh();
	viewer.data().clear();
	viewer.core.background_color = Eigen::Vector4f(1, 1, 1, 1);
	switched_mode = ((viewMode != prevMode) || (first_rendering));
	prevMode = viewMode;

	if ( (viewMode == ViewModeMesh) || (viewMode == ViewModeMeshWire) ) {
		render_mesh_and_wireframe(viewer);
	} else if (viewMode == ViewModeCreases) {
		render_crease_pattern(viewer);
	} else if (viewMode == ViewModeGauss) {
		render_gauss_map(viewer);
	} else if (viewMode == CreasesSVGReader) {
		render_crease_pattern_svg_reader(viewer);
	} else if (viewMode == ViewRulings) {
		render_mesh_and_wireframe(viewer);
		Eigen::MatrixXd VN; //igl::per_vertex_normals(Vren,Fren,VN);
    	igl::per_vertex_normals(dog->getV(),dog->getFTriangular(),VN);
		plot_vertex_based_rulings(viewer, dog->getV(), VN,
						dog->getQuadTopology(), new_rulings, rulings_length, rulings_mod, rulings_planar_eps);
	} else if (viewMode == ViewWallpaper) {
		render_wallpaper(viewer);
	}

	first_rendering = false;
}

void ModelViewer::render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = DC.getEditedSubmesh();
	if (switched_mode) viewer.core.align_camera_center(dog->getVrendering(), dog->getFrendering());
	if (render_curved_folding_properties) render_curved_folding_normals(viewer);
	//if ( state.dog.has_creases() && (DC.getEditedSubmeshI() <= -1) ) {
	render_dog_stitching_curves(viewer, state.dog, Eigen::RowVector3d(0, 0, 0));
	if (viewMode == ViewModeMeshWire) {
		if ( state.dog.has_creases() ) render_dog_wireframe(viewer);
		else render_wireframe(viewer, dog->getV(), dog->getQuadTopology()); 
	}
	if (render_pos_const) {
		render_positional_constraints(viewer);
		render_edge_points_constraints(viewer);
		DC.dogEditor->render_pairs();
	}
	render_mesh(viewer, dog->getVrendering(),dog->getFrendering());
}

void ModelViewer::render_wallpaper(igl::opengl::glfw::Viewer& viewer) {
	std::vector<Eigen::MatrixXd> Vlist; std::vector<Eigen::MatrixXi> Flist;
	const Dog* dog = DC.getEditedSubmesh(); Dog vertDog(*dog);
	//render_mesh(viewer,dog->getVrendering(),dog->getFrendering());
	auto left_curve = dog->left_bnd; auto right_curve = dog->right_bnd;
	auto lower_curve = dog->lower_bnd; auto upper_curve = dog->upper_bnd;
	Eigen::Matrix3d R; Eigen::Vector3d T;

	// add meshes to the right
	for (int j = 0; j < wallpaper_res; j++) {
		Vlist.push_back(vertDog.getVrendering()); Flist.push_back(vertDog.getFrendering());
		/*
		Dog nextDog(vertDog);
		for (int i = 0; i < wallpaper_res; i++) {
			PointsRigidAlignmentObjective::update_rigid_motion(nextDog.getV_vector(), left_curve, right_curve,R, T);
			Eigen::MatrixXd newV = (nextDog.getV() * R).rowwise() + T.transpose();
			nextDog.update_V(newV);
			Vlist.push_back(nextDog.getVrendering()); Flist.push_back(nextDog.getFrendering());
		}
		
		// add mesh up
		PointsRigidAlignmentObjective::update_rigid_motion(vertDog.getV_vector(), lower_curve, upper_curve,R, T);
		Eigen::MatrixXd newV = (vertDog.getV() * R).rowwise() + T.transpose();
		vertDog.update_V(newV);
		*/
	}
	

	Eigen::MatrixXd VWallpaper; Eigen::MatrixXi FWallpaper;
	igl::combine(Vlist,Flist, VWallpaper, FWallpaper);
	render_mesh(viewer,VWallpaper,FWallpaper);
}

void ModelViewer::render_crease_pattern(igl::opengl::glfw::Viewer& viewer) {
	if (switched_mode) viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
	int submesh_i = DC.getEditedSubmeshI();
	Dog flattenedDog(state.dog);
	if (switched_mode) viewer.core.align_camera_center(flattenedDog.getVrendering(), flattenedDog.getFrendering());
	if ( (submesh_i >= 0) && (submesh_i < state.dog.get_submesh_n()) ) {
		int submesh_v_min_i, submesh_v_max_i;
		flattenedDog.get_submesh_min_max_i(submesh_i, submesh_v_min_i, submesh_v_max_i, true);
		Eigen::MatrixXd& flattenedDogV = flattenedDog.getVMutable(); double z_offset = 1;
		for (int i = submesh_v_min_i; i <= submesh_v_max_i; i++) {
			flattenedDogV(i,2) += z_offset;
		}
		flattenedDog.update_Vren();
	}
	render_wireframe_boundary(viewer, flattenedDog.getV(), flattenedDog.getQuadTopology(), Eigen::RowVector3d(0.7, 0.7, 0.7));
	render_dog_stitching_constraints(viewer, flattenedDog, Eigen::RowVector3d(0, 0, 0.6));
	render_dog_stitching_curves(viewer, state.dog, Eigen::RowVector3d(0, 0, 0));
}

void ModelViewer::render_mesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& Vren, const Eigen::MatrixXi& Fren) {
	viewer.data().set_mesh(Vren, Fren);

	Eigen::Vector3d diffuse; diffuse << 135./255,206./255,250./255;
    Eigen::Vector3d ambient; /*ambient = 0.05*diffuse;*/ ambient<< 0.05,0.05,0.05;
    Eigen::Vector3d specular; specular << 0,0,0;
    //viewer.data.set_colors(diffuse);
    viewer.data().uniform_colors(ambient,diffuse,specular);

    Eigen::MatrixXd VN; igl::per_vertex_normals(Vren,Fren,VN);
  	viewer.data().set_normals(VN);
}

void ModelViewer::render_curved_folding_normals(igl::opengl::glfw::Viewer& viewer) {
	const DogEdgeStitching& eS = state.dog.getEdgeStitching();
	const Eigen::MatrixXd& V = state.dog.getV();
	for (int j = 0; j < eS.stitched_curves.size(); j++) {
		for (int i = 1; i < eS.stitched_curves[j].size()-1; i++) {
			EdgePoint eP =  eS.stitched_curves[j][i], eP_f = eS.stitched_curves[j][i+1], eP_b = eS.stitched_curves[j][i-1];
			Eigen::RowVector3d p0 = eP.getPositionInMesh(V), pf = eP_f.getPositionInMesh(V), pb = eP_b.getPositionInMesh(V);
			Eigen::RowVector3d t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
			Eigen::RowVector3d principal_n = ((pf-p0).normalized()+(pb-p0).normalized()).normalized();
			viewer.data().add_edges(p0, p0+principal_n, Eigen::RowVector3d(1,0,0));

			int v1,v2; state.dog.get_2_inner_vertices_from_edge(eP.edge,v1,v2);
			Eigen::RowVector3d t1 = V.row(v1)-p0, t2 = V.row(v2)-p0;
			// need to transpose when calling Eigen's cross product
			Eigen::RowVector3d n1 = -1*(t.cross(t1)).normalized();
			Eigen::RowVector3d n2 = -1*(t.cross(t2)).normalized();
			viewer.data().add_edges(p0, p0+n1, Eigen::RowVector3d(0,1,0));
			viewer.data().add_edges(p0, p0+n2, Eigen::RowVector3d(0,0,1));

			//std::cout << "cos angle t,n1 = " << t.dot(n1) << std::endl;

			auto angle1 = acos(principal_n.dot(n1))*180.0 / M_PI, angle2 = acos(principal_n.dot(n2))*180.0 / M_PI;
			//std::cout << "angles at i = " << i << " are " << angle1 << ", " << angle2 << " with abs(diff) = " << abs(angle1-angle2) << std::endl;
		}
	}
	cout << endl << endl;
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
	const Dog* dog = DC.getEditedSubmesh();
	Eigen::VectorXd x(dog->getV_vector());
	Eigen::VectorXi b; Eigen::VectorXd bc; DC.get_positional_constraints(b,bc);
	Eigen::VectorXd constrained_pts_coords_vec; igl::slice(x,b,1, constrained_pts_coords_vec);

	int pts_num = b.size()/3;
	Eigen::MatrixXd E1(pts_num,3),E2(pts_num,3);
	for (int i = 0; i < pts_num; i++) {
		E1.row(i) << constrained_pts_coords_vec(i),constrained_pts_coords_vec(pts_num+i),constrained_pts_coords_vec(2*pts_num+i);
		E2.row(i) << bc(i),bc(pts_num+i),bc(2*pts_num+i);
	}
	
	viewer.data().add_edges(E1,E2,Eigen::RowVector3d(1.,0,0));
	//viewer.data().add_points(E1,Eigen::RowVector3d(1.,0,0));
	DC.dogEditor->render_positional_constraints();
	for (int d_i = 0; d_i < DC.dihedral_constrained.size(); d_i++) {
		auto ep = DC.dihedral_constrained[d_i];
		viewer.data().add_points(ep.getPositionInMesh(DC.getEditedSubmesh()->getV()),Eigen::RowVector3d(1.,1,0));
	}
}

void ModelViewer::render_edge_points_constraints(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = DC.getEditedSubmesh();
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
	DC.get_edge_point_constraints(edgePoints, edgeCoords);
	if (!edgePoints.size()) return;
	Eigen::MatrixXd currentCoords = EdgePoint::getPositionInMesh(edgePoints, dog->getV());
	viewer.data().add_edges(currentCoords,edgeCoords,Eigen::RowVector3d(1.,0,0));
	//viewer.data().add_points(currentCoords,Eigen::RowVector3d(1.,0,0));
	//viewer.data().add_points(edgeCoords,Eigen::RowVector3d(1.,0,0));
}

void ModelViewer::render_gauss_map(igl::opengl::glfw::Viewer& viewer) {
  const Dog* dog = DC.getEditedSubmesh();
  viewer.data().set_mesh(sphereV, sphereF);
  Eigen::Vector3d diffuse; diffuse << 0.98,0.98,0.98;
  Eigen::Vector3d ambient; ambient << 0,0,0;//0.05*diffuse;
  Eigen::Vector3d specular; specular << 0.05*diffuse;
  viewer.data().uniform_colors(ambient,diffuse,specular);
  //viewer.core.shininess = 0;
  
  if (switched_mode) viewer.core.align_camera_center(dog->getVrendering(), dog->getFrendering());
  //viewer.core.align_camera_center(sphereV, sphereF);
  //viewer.core.show_lines = false;

  // TODO support curved folds by looking at normal map of each one separately
  Eigen::MatrixXd VN; igl::per_vertex_normals(dog->getVrendering(),dog->getFrendering(),VN);
  //viewer.data.set_normals(VN);
  render_wireframe(viewer,VN,dog->getQuadTopology(), false);
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

void ModelViewer::render_dog_wireframe(igl::opengl::glfw::Viewer& viewer) {
	const std::vector<std::pair<int,int>>& wire_edges_i =  DC.getEditedSubmesh()->get_rendered_wireframe_edges();
	const Eigen::MatrixXd& V_ren =  DC.getEditedSubmesh()->getVrendering();
	int e_num = wire_edges_i.size();
	Eigen::MatrixXd E1(e_num,3), E2(e_num,3);
	for (int i = 0; i < e_num; i++) {
		E1.row(i) = V_ren.row(wire_edges_i[i].first);
		E2.row(i) = V_ren.row(wire_edges_i[i].second);
	}
	viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0,0,0));
}