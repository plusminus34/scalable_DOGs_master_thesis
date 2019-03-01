#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "Lasso.h"

#include "../Dog/Dog.h"
//#include "Editor.h"

class DogEditor {
public:
	enum EditMode { SELECT_POSITIONAL, TRANSLATE, VERTEX_PAIRS, EDGES_ANGLE, DIHEDRAL_ANGLE, NONE};//, ROTATE, CUT, GLUE1, GLUE2, NONE };
	enum SelectMode {VertexPicker, EdgePointPicker, CurvePicker};

	DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, EditMode& edit_mode, SelectMode& select_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	~DogEditor();

	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up();

	void reset_new_constraint();

	void render_pairs() const {render_paired_constraints();render_selected_pairs();}

	void render_positional_constraints() const;
	void render_paired_constraints() const;
	void render_selected_pairs() const;
	void clearHandles();

	int pair_vertex_1 = -1; int pair_vertex_2 = -1; bool next_pair_first = true;
	int edge_angle_v1 = -1; int edge_angle_center = -1; int edge_angle_v2 = -1; int edges_angle_pick_idx = 0;

	EdgePoint picked_edge;
private:
	Dog& dog;
	EditMode& edit_mode; SelectMode& select_mode;
	bool& has_new_constraints;
	Eigen::VectorXi& b; Eigen::VectorXd& bc; std::vector<std::pair<int,int>>& paired_vertices;
	std::vector<EdgePoint>& edgePoints; Eigen::MatrixXd& edgeCoords;
	
	// TODO: Use V_ren and F_ren with the editor. Convert inner points to the correct V index and edge points to edge point constraints
	int pick_vertex();
	int pick_edge(EdgePoint& edgePoint);

	void select_positional_mouse_down();
	void translate_vertex_edit_mouse_down();
	void vertex_pairs_edit_mouse_down();
	void edges_angle_edit_mouse_down();
	void dihedral_angle_edit_mouse_down();

	void onNewHandleID();
	void compute_handle_centroids();
	void get_new_handle_locations(Eigen::Vector3f translation);

	void reset_new_pair_constraint();
	void reset_new_edge_angle_constraint();

	igl::opengl::glfw::Viewer& viewer;
	const Eigen::MatrixXd &V_ren; const Eigen::MatrixXd &V; Eigen::MatrixXd oldV;
	Eigen::MatrixXi F; // The triangular faces of the dog (not dogFrendering)
	Lasso lasso;

	//Eigen::Vector3f translation;

	//updated positions of handle vertices, #HV x3
	Eigen::MatrixXd handle_vertex_positions;//(0,3);
	//list of all vertices belonging to handles, #HV x1
	Eigen::VectorXi handle_vertices;//(0,1);
	//for saving constrained vertices
	//vertex-to-handle index, #V x1 (-1 if vertex is free)
	Eigen::VectorXi handle_id;//(0,1);
	int current_handle = -1;

	//index of handle being moved
	int moving_handle = -1;

	//centroids of handle regions, #H x1
	Eigen::MatrixXd handle_centroids;//(0,3);
	int down_mouse_x = -1, down_mouse_y = -1;

	int spline_pt_number = 4;

	bool action_started = false;
};