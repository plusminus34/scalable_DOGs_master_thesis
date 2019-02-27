#include "DogEditor.h"

#include <igl/project.h>
#include <igl/quat_conjugate.h>
#include <igl/quat_mult.h>
#include <igl/slice.h>
#include <igl/unproject.h>

using namespace std;

std::vector<int> get_second_dog_row(Dog& dog);

DogEditor::DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, MouseMode& mouse_mode, SelectMode& select_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : 
		dog(dog), viewer(viewer),
		has_new_constraints(has_new_constraints),b(b),bc(bc),
		paired_vertices(paired_vertices),edgePoints(edgePoints), edgeCoords(edgeCoords),
		V(dog.getV()), F(dog.getFTriangular()),
		lasso(viewer,V,F), mouse_mode(mouse_mode), select_mode(select_mode) {
	
	handle_id.setConstant(V.rows(), 1, -1);
	oldV = V;
}

/*
std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}*/

bool DogEditor::callback_mouse_down() {
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;
	
	if (mouse_mode == SELECT) {
		if (select_mode == CurvePicker) {
			if (lasso.strokeAddCurve(viewer.current_mouse_x, viewer.current_mouse_y) >=0) {
				action_started = true;
			}
		} else if (select_mode == VertexPicker) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if (vi >=0) {
				int index = handle_id.maxCoeff()+1;
				if (handle_id[vi] == -1) handle_id[vi] = index;
				current_handle = index;
				
				onNewHandleID();
			}
		} else if (select_mode == PairPicker) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if (vi >=0) {
				if (next_pair_first) {
					pair_vertex_1 = vi; next_pair_first = false;
				} else {
					pair_vertex_2 = vi; next_pair_first = true;
				}
			}
		}
	} else if (mouse_mode == TRANSLATE) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if(vi>=0 && handle_id[vi]>=0)  {//if a region was found, mark it for translation/rotation {
				moving_handle = handle_id[vi];
				current_handle = moving_handle;
				oldV = V;//.copy();
				//get_new_handle_locations();
				//compute_grad_constraints();
				action_started = true;
			}
	} else if (mouse_mode == APPLY) {
		if (select_mode == PairPicker) {
			int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
			if ( (pair_vertex_1!=-1) && (pair_vertex_2!=-1) ){
				if ((vi == pair_vertex_1) || (vi == pair_vertex_2)) {
					int vnum = V.rows();
					for (int i = 0; i < 3; i++) {
						paired_vertices.push_back(std::pair<int,int>(i*vnum+pair_vertex_1,i*vnum+pair_vertex_2));	
					}
					
					pair_vertex_1 = pair_vertex_2 = -1;
					has_new_constraints = true;
				}
			}
		}
	}
	return action_started;
}
bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	if (!action_started)
		return false;
	
	if (mouse_mode == SELECT) {
		if (select_mode == CurvePicker) {
			lasso.strokeAddCurve(mouse_x, mouse_y);
			return true;
		}
	}
	
	if (mouse_mode == TRANSLATE) {
		Eigen::Vector3f translation = computeTranslation(mouse_x, down_mouse_x, mouse_y, down_mouse_y, handle_centroids.row(moving_handle));
		get_new_handle_locations(translation);
		down_mouse_x = mouse_x;
		down_mouse_y = mouse_y;
		return true;
	}
	return false;
}

bool DogEditor::callback_mouse_up() {
	if (!action_started)
		return false;
	action_started = false;
	
	if (mouse_mode == TRANSLATE) {
		Eigen::Vector3f translation; translation.setZero();
		moving_handle = -1;
		oldV = V;
		
		lasso.reinit();
		compute_handle_centroids();
		get_new_handle_locations(translation);
		
		return true;
	}

	if (mouse_mode == SELECT) {
		if (select_mode == CurvePicker) {
			lasso.strokeFinishCurve(spline_pt_number);
			//lasso->strokeFinishCurve(D.const_edges, D.edge_coordinates);
			//lasso->strokeFacesSqr.clear();
		}
		return true;
	}
	return false;
};

void DogEditor::onNewHandleID() {
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
	b.resize(3*handle_vertices.rows()); bc.resize(b.rows());
	for (int i = 0; i < const_v_num; i++) {
		b(i) = handle_vertices(i); bc(i) = handle_vertex_positions(i,0);
		b(const_v_num+i) = handle_vertices(i) + v_num; bc(const_v_num+i) = handle_vertex_positions(i,1);
		b(2*const_v_num+i) = handle_vertices(i) + 2*v_num; bc(2*const_v_num+i) = handle_vertex_positions(i,2);
	}
	//solver.setConst(D.b);
	//Direction direction; direction.x = direction.y = direction.z = 0;
	//handle_directions.push_back(direction);
	has_new_constraints = true;
}

void DogEditor::get_new_handle_locations(Eigen::Vector3f translation) {
	int count = 0;
	for (long vi = 0; vi<V.rows(); ++vi)
		if(handle_id[vi] >=0) {
			//Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
			Eigen::RowVector3f goalPosition = oldV.row(vi).cast<float>();
			//Eigen::RowVector3f goalPosition = oldV.row(vi).cast<float>();
			if (handle_id[vi] == moving_handle){
				if( mouse_mode == TRANSLATE) goalPosition+=translation;
			}
			handle_vertex_positions.row(count++) = goalPosition.cast<double>();
			oldV.row(vi) = goalPosition.cast<double>();;
		}
	const int const_v_num = handle_vertices.rows(); const int v_num = V.rows();
	bc.resize(b.rows());
	for (int i = 0; i < const_v_num; i++) {
		bc(i) = handle_vertex_positions(i,0);
		bc(const_v_num+i) = handle_vertex_positions(i,1);
		bc(2*const_v_num+i) = handle_vertex_positions(i,2);
	}
}

void DogEditor::compute_handle_centroids() {
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
Eigen::Vector3f DogEditor::computeTranslation(int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D) {
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
Eigen::Vector4f DogEditor::computeRotation(int mouse_x, int from_x,int mouse_y, int from_y, Eigen::RowVector3d pt3D) {
	
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

void DogEditor::render_positional_constraints() const {
	Eigen::MatrixXd const_v;
	Eigen::MatrixXd handle_colors(handle_vertex_positions.rows(),3);
    for (int i = 0; i < handle_colors.rows(); i++) {
      if (handle_id[handle_vertices[i]] == moving_handle) {
        handle_colors.row(i) = Eigen::RowVector3d(220./255,0./255,102./255);
      } else {
        handle_colors.row(i) = Eigen::RowVector3d(30./255,80./255,255./255);
      }
    }
    igl::slice(V, handle_vertices, 1, const_v);
    viewer.data().add_points(const_v, handle_colors);
}
void DogEditor::render_paired_constraints() const {
	Eigen::MatrixXd E1(paired_vertices.size()/3,3),E2(paired_vertices.size()/3,3);
	// The pair values are flattened
	for (int i = 0; i < paired_vertices.size(); i+=3) {
		E1.row(i/3) = V.row(paired_vertices[i].first);
		E2.row(i/3) = V.row(paired_vertices[i].second);
	}
	viewer.data().add_edges(E1, E2, Eigen::RowVector3d(128./255,128./255,128./255));
	viewer.data().add_points(E1, Eigen::RowVector3d(128./255,128./255,128./255));
	viewer.data().add_points(E2, Eigen::RowVector3d(128./255,128./255,128./255));
}
void DogEditor::render_selected_pairs() const {
	Eigen::RowVector3d active_pair_color(0./255,160./255,0./255);
	if (pair_vertex_1 != -1) {
		Eigen::RowVector3d pt1 = V.row(pair_vertex_1);
		viewer.data().add_points(pt1, active_pair_color);	
	}
	if (pair_vertex_2 != -1) {
		Eigen::RowVector3d pt2 = V.row(pair_vertex_2);
		viewer.data().add_points(pt2, active_pair_color);	
	}	
}

void DogEditor::clearHandles() {
	handle_id.setConstant(V.rows(),1,-1);
	handle_vertex_positions.setZero(0,3);
	handle_vertices.resize(0);handle_vertices.setZero(0);

	moving_handle = -1;
	current_handle = -1;
	pair_vertex_1 = pair_vertex_2 = -1;
	paired_vertices.clear();
}