#include "DogEditor.h"

#include "EditingUtils.h"
#include <igl/slice.h>

using namespace std;

std::vector<int> get_second_dog_row(Dog& dog);

DogEditor::DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, EditMode& edit_mode, SelectMode& select_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : 
		dog(dog), viewer(viewer),
		has_new_constraints(has_new_constraints),b(b),bc(bc),
		paired_vertices(paired_vertices),edgePoints(edgePoints), edgeCoords(edgeCoords),
		V(dog.getV()), F(dog.getFTriangular()),
		lasso(viewer,V,F), edit_mode(edit_mode), select_mode(select_mode) {
	
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
	
	cout << "mouse down, and edit_mode = " << edit_mode <<  endl;
	if (edit_mode == SELECT_POSITIONAL) {
		select_positional_mouse_down();
	} else if (edit_mode == TRANSLATE) {
		translate_vertex_edit_mouse_down();
	} else if (edit_mode == VERTEX_PAIRS) {
		vertex_pairs_edit_mouse_down();
	} else if (edit_mode == EDGES_ANGLE) {
		edges_angle_edit_mouse_down();
	}
	return action_started;
}

bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	if (!action_started)
		return false;
	
	if (edit_mode == SELECT_POSITIONAL) {
		if (select_mode == CurvePicker) {
			lasso.strokeAddCurve(mouse_x, mouse_y);
			return true;
		}
	} else if (edit_mode == TRANSLATE) {
		Eigen::Vector3f translation = computeTranslation(viewer, mouse_x, down_mouse_x, mouse_y, down_mouse_y, handle_centroids.row(moving_handle));
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
	
	if (edit_mode == TRANSLATE) {
		Eigen::Vector3f translation; translation.setZero();
		moving_handle = -1;
		oldV = V;
		
		lasso.reinit();
		compute_handle_centroids();
		get_new_handle_locations(translation);
		
		return true;
	}

	if (edit_mode == SELECT_POSITIONAL) {
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

	// b contains normal constraints
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
				if( edit_mode == TRANSLATE) goalPosition+=translation;
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

void DogEditor::apply_new_constraint() {
	if (edit_mode == VERTEX_PAIRS) {
		if ( (pair_vertex_1!= -1) && (pair_vertex_2!= -1) ) {
			int vnum = V.rows();
			for (int i = 0; i < 3; i++) {
				paired_vertices.push_back(std::pair<int,int>(i*vnum+pair_vertex_1,i*vnum+pair_vertex_2));	
			}
		
			pair_vertex_1 = pair_vertex_2 = -1;
			has_new_constraints = true;	
		}
	}
}

void DogEditor::cancel_new_constraint() {
	if (edit_mode == VERTEX_PAIRS) {
		
	} else if (edit_mode == EDGES_ANGLE) {
		edge_angle_v1 = edge_angle_v2 = edge_angle_center = -1; edges_angle_pick_idx = 0;
	}
}

void DogEditor::cancel_new_pair_constraint() {
	pair_vertex_1 = pair_vertex_2 = -1; next_pair_first = true;
}

void DogEditor::cancel_new_edge_angle_constraint() {
	edge_angle_v1 = -1; edge_angle_center = -1; edge_angle_v2 = -1; edges_angle_pick_idx = 0;
}

void DogEditor::clearHandles() {
	handle_id.setConstant(V.rows(),1,-1);
	handle_vertex_positions.setZero(0,3);
	handle_vertices.resize(0); handle_vertices.setZero(0);

	moving_handle = -1;
	current_handle = -1;
	cancel_new_pair_constraint(); paired_vertices.clear();
	cancel_new_edge_angle_constraint();
}

void DogEditor::select_positional_mouse_down() {
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
	}
}

void DogEditor::translate_vertex_edit_mouse_down() {
	int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
	if(vi>=0 && handle_id[vi]>=0)  {//if a region was found, mark it for translation/rotation {
		moving_handle = handle_id[vi];
		current_handle = moving_handle;
		oldV = V;//.copy();
		//get_new_handle_locations();
		//compute_grad_constraints();
		action_started = true;
	}
}

void DogEditor::vertex_pairs_edit_mouse_down() {
	if (select_mode != VertexPicker) return;
	int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
	if (vi >=0) {
		if (next_pair_first) {
			pair_vertex_1 = vi; next_pair_first = false;
		} else {
			pair_vertex_2 = vi; next_pair_first = true;
		}
	}
}

void DogEditor::edges_angle_edit_mouse_down() {
	if (select_mode != VertexPicker) return;
	int vi = lasso.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
	if (vi >=0) {
		if (edges_angle_pick_idx == 0) {
			edge_angle_v1 = vi; edges_angle_pick_idx +=1;
		} else if (edges_angle_pick_idx == 1) {
			edge_angle_center = vi; edges_angle_pick_idx +=1;
		} else {
			edge_angle_v2 = vi; edges_angle_pick_idx +=1;
		}
	}
	edges_angle_pick_idx = edges_angle_pick_idx % 3;
}