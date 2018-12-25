#include "Editor.h"

#include <igl/project.h>
#include <igl/quat_conjugate.h>
#include <igl/quat_mult.h>
#include <igl/slice.h>
#include <igl/unproject.h>

using namespace std;

Editor::Editor(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F_tri) : 
				viewer(viewer), V(V), F(F_tri), lasso(viewer,V,F), translation(0,0,0) {}

bool Editor::callback_mouse_down(int button, int modifier) {
	
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;
	
	if (mouse_mode == SELECT) {
		if (select_mode == CurvePicker) {
			// TODO
		} else if (select_mode == VertexPicker) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if (vi >=0) {
				selected_v.resize(1,1);
				selected_v[0] = vi;
				applySelection();
			}
		} else if (select_mode == PathPicker) {
			// TODO
		}
	} else if (mouse_mode == TRANSLATE) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if(vi>=0 && handle_id[vi]>=0)  {//if a region was found, mark it for translation/rotation {
				moving_handle = handle_id[vi];
				current_handle = moving_handle;
				get_new_handle_locations();
				//compute_grad_constraints();
				deforming = true;
			}
	}
	return deforming;
}

void Editor::applySelection() {
	int index = handle_id.maxCoeff()+1;
	for (int i =0; i < selected_v.rows(); ++i)
	{
		const int selected_vertex = selected_v[i];
		cout << "Selected V  " << selected_vertex << " at location " << V.row(selected_vertex) << endl;
		if (handle_id[selected_vertex] == -1)
			handle_id[selected_vertex] = index;
	}
	selected_v.resize(0,1);
	current_handle = index;
	
	onNewHandleID();
}

void Editor::onNewHandleID() {
	//store handle vertices too
	int numFree = (handle_id.array()==-1).cast<int>().sum();
	int num_handle_vertices = V.rows() - numFree;
	handle_vertices.setZero(num_handle_vertices);
	handle_vertex_positions.setZero(num_handle_vertices,3);
	
	int count = 0;
	for (long vi = 0; vi<V.rows(); ++vi)
		if(handle_id[vi] >=0)
			handle_vertices[count++] = vi;
	
	compute_handle_centroids();

	// update handle_vertices_fixed_pos
	igl::slice(V, handle_vertices, 1, handle_vertex_positions);
	// update b and bc (vector form of handle_vertices and handle_vertex_positions)
	const int const_v_num = handle_vertices.rows(); const int v_num = V.rows();

	// b contains normal constraints + z-up/down constraints for folding
	b.resize(3*handle_vertices.rows());
	for (int i = 0; i < const_v_num; i++) {
		b(3*i) = handle_vertices(i); bc(3*i) = handle_vertex_positions(i,0);
		b(3*i+1) = handle_vertices(i) + v_num; bc(3*i+1) = handle_vertex_positions(i,1);
		b(3*i+2) = handle_vertices(i) + 2*v_num; bc(3*i+2) = handle_vertex_positions(i,2);
	}
	//solver.setConst(D.b);
	//Direction direction; direction.x = direction.y = direction.z = 0;
	//handle_directions.push_back(direction);
}

void Editor::get_new_handle_locations() {
	int count = 0;
	for (long vi = 0; vi<V.rows(); ++vi)
		if(handle_id[vi] >=0) {
			Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
			//Eigen::RowVector3f goalPosition = oldV.row(vi).cast<float>();
			if (handle_id[vi] == moving_handle){
				if( mouse_mode == TRANSLATE) goalPosition+=translation;
			}
			handle_vertex_positions.row(count++) = goalPosition.cast<double>();
			//D.oldV.row(vi) = goalPosition.cast<double>();;
		}
	const int const_v_num = handle_vertices.rows(); const int v_num = V.rows();
	bc.resize(b.rows());
	for (int i = 0; i < const_v_num; i++) {
		bc(3*i) = handle_vertex_positions(i,0);
		bc(3*i+1) = handle_vertex_positions(i,1);
		bc(3*i+2) = handle_vertex_positions(i,2);
	}
}

void Editor::compute_handle_centroids() {
	//compute centroids of handles
	int num_handles = handle_id.maxCoeff()+1;
	handle_centroids.setZero(num_handles,3);
	
	Eigen::VectorXi num; num.setZero(num_handles,1);
	for (long vi = 0; vi<V.rows(); ++vi)
	{
		int r = handle_id[vi];
		if ( r!= -1)
		{
			handle_centroids.row(r) += V.row(vi);
			num[r]++;
		}
	}
	
	for (long i = 0; i<num_handles; ++i)
		handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];
	
}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f Editor::computeTranslation(int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D) {
	//Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::Matrix4f modelview = viewer.core.view;// * viewer.core.model;
	//project the given point (typically the handle centroid) to get a screen space depth
	Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
																			modelview,
																			viewer.core.proj,
																			viewer.core.viewport);
	float depth = proj[2];
	
	double x, y;
	Eigen::Vector3f pos1, pos0;
	
	//unproject from- and to- points
	x = mouse_x;
	y = viewer.core.viewport(3) - mouse_y;
	pos1 = igl::unproject(Eigen::Vector3f(x,y,depth),
												modelview,
												viewer.core.proj,
												viewer.core.viewport);
	
	
	x = from_x;
	y = viewer.core.viewport(3) - from_y;
	pos0 = igl::unproject(Eigen::Vector3f(x,y,depth),
												modelview,
												viewer.core.proj,
												viewer.core.viewport);
	
	//translation is the vector connecting the two
	Eigen::Vector3f translation;
	translation = pos1 - pos0;
	
	return translation;
}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f Editor::computeRotation(int mouse_x, int from_x,int mouse_y, int from_y, Eigen::RowVector3d pt3D) {
	
	Eigen::Vector4f rotation;
	rotation.setZero();
	rotation[3] = 1.;
	
	Eigen::Matrix4f modelview = viewer.core.view;// * viewer.core.model;
	
	//initialize a trackball around the handle that is being rotated
	//the trackball has (approximately) width w and height h
	double w = viewer.core.viewport[2]/8;
	double h = viewer.core.viewport[3]/8;
	
	//the mouse motion has to be expressed with respect to its center of mass
	//(i.e. it should approximately fall inside the region of the trackball)
	
	//project the given point on the handle(centroid)
	Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
																			modelview,
																			viewer.core.proj,
																			viewer.core.viewport);
	proj[1] = viewer.core.viewport[3] - proj[1];
	
	//express the mouse points w.r.t the centroid
	from_x -= proj[0]; mouse_x -= proj[0];
	from_y -= proj[1]; mouse_y -= proj[1];
	
	//shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
	from_x += w/2; mouse_x += w/2;
	from_y += h/2; mouse_y += h/2;
	
	//get rotation from trackball
	//Eigen::Vector4f drot = viewer.core.trackball_angle;
	Eigen::Quaternionf drot_t = viewer.core.trackball_angle;
	Eigen::Vector4f drot = drot_t.coeffs();
	Eigen::Vector4f drot_conj;
	igl::quat_conjugate(drot.data(), drot_conj.data());
	igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());
	
	//account for the modelview rotation: prerotate by modelview (place model back to the original
	//unrotated frame), postrotate by inverse modelview
	Eigen::Vector4f out;
	igl::quat_mult(rotation.data(), drot.data(), out.data());
	igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
	return rotation;
}