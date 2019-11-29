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

	igl::readOBJ("../../data/sphere.obj",sphereV,sphereF);
	center_and_scale_gauss_sphere(sphereV,sphereF);
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = DC.getEditedSubmesh();
	clear_edges_and_points(viewer);
	if ((viewMode != prevMode) || (first_rendering)) switched_mode = true;
	prevMode = viewMode;
	if (first_rendering || switched_mode)  {
		viewer.data().clear();
		viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 1);
	}

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
	} else if (viewMode == ViewCurves){
		render_curves(viewer);
	}

	first_rendering = false;
}

void ModelViewer::clear_edges_and_points(igl::opengl::glfw::Viewer& viewer) {
	viewer.data().set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
	viewer.data().set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
}

void ModelViewer::render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = DC.getEditedSubmesh();
	if (switched_mode) viewer.core().align_camera_center(dog->getVrendering(), dog->getFrendering());
	if (render_curved_folding_properties) render_curved_osculating_planes(viewer);
	//if ( state.dog.has_creases() && (DC.getEditedSubmeshI() <= -1) ) {
	if (show_curves) render_dog_stitching_curves(viewer, state.dog, Eigen::RowVector3d(0, 0, 0));
	if (viewMode == ViewModeMeshWire) {
		if ( state.dog.has_creases() && culled_view) render_dog_wireframe(viewer);
		else render_wireframe(viewer, dog->getV(), dog->getQuadTopology());
	}
	if (show_conversion) render_conversion(viewer);
	render_positional_constraints(viewer);
	DC.dogEditor->render_pairs();
	if (render_pos_const) {
		render_edge_points_constraints(viewer);
	}
	if (culled_view) render_mesh(viewer, dog->getVrendering(),dog->getFrendering()); else render_mesh(viewer, dog->getV(),dog->getFTriangular());
}

void ModelViewer::render_crease_pattern(igl::opengl::glfw::Viewer& viewer) {
	if (switched_mode) viewer.core().align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
	int submesh_i = DC.getEditedSubmeshI();
	Dog flattenedDog(state.dog);
	if (switched_mode) viewer.core().align_camera_center(flattenedDog.getVrendering(), flattenedDog.getFrendering());
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
	if (show_curves) render_dog_stitching_curves(viewer, state.dog, Eigen::RowVector3d(0, 0, 0));
}

void ModelViewer::render_mesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& Vren, const Eigen::MatrixXi& Fren) {
	if (first_rendering || switched_mode) {
		viewer.data().set_mesh(Vren, Fren);
		//Eigen::Vector3d diffuse; diffuse << 135./255,206./255,250./255;
		Eigen::Vector3d diffuse; diffuse << 0,0,0;
	    Eigen::Vector3d ambient; /*ambient = 0.05*diffuse;*/ ambient<< 210.0/255,237.0/255,1.0;
	    Eigen::Vector3d specular; specular << 0,0,0;
	    //viewer.data.set_colors(diffuse);
	    viewer.data().uniform_colors(ambient,diffuse,specular);
	}
	else {
		 viewer.data().set_vertices(Vren);
    	 viewer.data().compute_normals();
	}

    //Eigen::MatrixXd VN; igl::per_vertex_normals(Vren,Fren,VN);
  	//viewer.data().set_normals(VN);
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

void ModelViewer::render_curved_osculating_planes(igl::opengl::glfw::Viewer& viewer) {
	const DogEdgeStitching& eS = state.dog.getEdgeStitching();
	const Eigen::MatrixXd& V = state.dog.getV();
	for (int j = 0; j < eS.stitched_curves.size(); j++) {
		for (int i = 1; i < eS.stitched_curves[j].size()-1; i+=2) {
			EdgePoint eP =  eS.stitched_curves[j][i], eP_f = eS.stitched_curves[j][i+1], eP_b = eS.stitched_curves[j][i-1];
			Eigen::RowVector3d p0 = eP.getPositionInMesh(V), pf = eP_f.getPositionInMesh(V), pb = eP_b.getPositionInMesh(V);

			Eigen::RowVector3d osculating_plane_n = (p0-pf).cross(p0-pb).normalized();
			Eigen::RowVector3d t1 = (p0-pf).cross(osculating_plane_n).normalized();
			Eigen::RowVector3d t2 = (t1).cross(osculating_plane_n).normalized();

			double half_l = 0.4; Eigen::RowVector3d plane_color(155./255,198./255,161./255);
			Eigen::RowVector3d corner_1 = p0+half_l*t1+half_l*t2, corner_2 = p0+half_l*t1-half_l*t2, corner_3 = p0-half_l*t1-half_l*t2, corner_4 = p0-half_l*t1+half_l*t2;
			viewer.data().add_edges(corner_1, corner_2, plane_color);
			viewer.data().add_edges(corner_2, corner_3, plane_color);
			viewer.data().add_edges(corner_3, corner_4, plane_color);
			viewer.data().add_edges(corner_4, corner_1, plane_color);
		}
	}
}

void ModelViewer::render_crease_pattern_svg_reader(igl::opengl::glfw::Viewer& viewer) {
	if (state.dog.has_creases()) {
		//viewer.data().set_mesh(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
		viewer.data().set_face_based(true);
		viewer.data().set_colors(state.creasesVisualization.faceColors);
		viewer.data().add_edges(state.creasesVisualization.edge_pts1,state.creasesVisualization.edge_pts2,Eigen::RowVector3d(0,0,0));
		viewer.data().add_edges(state.creasesVisualization.meshE1,state.creasesVisualization.meshE2,Eigen::RowVector3d(0,0,0));
		//viewer.core().align_camera_center(state.creasesVisualization.V_arr, state.creasesVisualization.F_arr);
		viewer.core().align_camera_center(state.creasesVisualization.edge_pts1);
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

	if (render_pos_const) viewer.data().add_edges(E1,E2,Eigen::RowVector3d(1.,0,0));
	//viewer.data().add_points(E1,Eigen::RowVector3d(1.,0,0));
	DC.dogEditor->render_positional_constraints();
	for (int d_i = 0; d_i < DC.dihedral_constrained.size(); d_i++) {
		auto ep = DC.dihedral_constrained[d_i];
		viewer.data().add_points(ep.getPositionInMesh(DC.getEditedSubmesh()->getV()),Eigen::RowVector3d(155./255,198./255,161./255));
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
  //viewer.core().shininess = 0;

  if (switched_mode) viewer.core().align_camera_center(dog->getVrendering(), dog->getFrendering());
  //viewer.core().align_camera_center(sphereV, sphereF);
  //viewer.core().show_lines = false;

  // TODO support curved folds by looking at normal map of each one separately
  Eigen::MatrixXd VN; igl::per_vertex_normals(dog->getVrendering(),dog->getFrendering(),VN);
  //viewer.data.set_normals(VN);
  render_wireframe(viewer,VN,dog->getQuadTopology(), false);
  if (switched_mode) viewer.core().align_camera_center(sphereV, sphereF);
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

void ModelViewer::render_conversion(igl::opengl::glfw::Viewer& viewer) {
	Eigen::MatrixXd V = state.dog.getV();
	Eigen::MatrixXd colors(V.rows(), 3);
	for(int i=0; i<V.rows(); ++i){
		int c = state.conversion.fine_to_coarse(i);
		if(c > -1){
			colors.row(i) << 0.0, 0.7, 0.3;
		} else if (c == -2) {
			colors.row(i) << 0.0, 0.0, 0.4;
		} else {
			colors.row(i) << 0.3, 0.0, 0.0;
		}
	}
	viewer.data().add_points(V, colors);
}

void ModelViewer::render_curves(igl::opengl::glfw::Viewer& viewer) {
	// get coarse curve
	const Eigen::MatrixXd& coarse_V = state.coarse_dog.getV() *2;//coarse scale
	const DogEdgeStitching& coarse_eS = state.coarse_dog.getEdgeStitching();
	int num_curves = coarse_eS.stitched_curves.size();
	int num_edges = -num_curves;
	for(int i=0; i<num_curves; ++i) num_edges += coarse_eS.stitched_curves[i].size();
	Eigen::MatrixXd coarse_E1(num_edges, 3), coarse_E2(num_edges, 3);
	int ie = 0;
	for (int j = 0; j < num_curves; ++j) {
		Eigen::RowVector3d p0 = coarse_eS.stitched_curves[j][0].getPositionInMesh(coarse_V);
		for(int k=1; k<coarse_eS.stitched_curves[j].size(); ++k){
			Eigen::RowVector3d p1 = coarse_eS.stitched_curves[j][k].getPositionInMesh(coarse_V);
			coarse_E1.row(ie) = p0;
			coarse_E2.row(ie) = p1;
			++ie;
			p0 = p1;
		}
	}
	// get fine curve (aka the one you see normally)
	const Eigen::MatrixXd& V = state.dog.getV();
	const DogEdgeStitching& eS = state.dog.getEdgeStitching();
	num_edges = -num_curves;
	for(int i=0; i<num_curves; ++i) num_edges += eS.stitched_curves[i].size();
	Eigen::MatrixXd E1(num_edges, 3), E2(num_edges, 3);
	ie = 0;
	for (int j = 0; j < num_curves; ++j) {
		Eigen::RowVector3d p0 = eS.stitched_curves[j][0].getPositionInMesh(V);
		for(int k=1; k<eS.stitched_curves[j].size(); ++k){
			Eigen::RowVector3d p1 = eS.stitched_curves[j][k].getPositionInMesh(V);
			E1.row(ie) = p0;
			E2.row(ie) = p1;
			++ie;
			p0 = p1;
		}
	}
	// get interpolated curve
	Eigen::MatrixXd approx_E1(num_edges, 3), approx_E2(num_edges, 3);
	ie = 0;
	for (int j = 0; j < num_curves; ++j) {
		Eigen::MatrixXd interpolated_points = state.conversion.getInterpolatedCurveCoords(state.dog, state.coarse_dog, j);
		for(int k=1; k<interpolated_points.rows(); ++k){
			approx_E1.row(ie) = interpolated_points.row(k-1) *2;//coarse scale
			approx_E2.row(ie) = interpolated_points.row(k) *2;//coarse scale
			++ie;
		}
	}
	// figure out what to display: Everything
	viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0,0,0));
	viewer.data().add_edges(coarse_E1, coarse_E2, Eigen::RowVector3d(0.75,0,0));
	viewer.data().add_edges(approx_E1, approx_E2, Eigen::RowVector3d(0.25,0.75,0.25));
}
