#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "Lasso.h"

#include "../Dog/Dog.h"
//#include "Editor.h"

class DogEditor {
public:
	enum MouseMode { SELECT, TRANSLATE, APPLY, NONE};//, ROTATE, CUT, GLUE1, GLUE2, NONE };
	enum SelectMode {VertexPicker, PairPicker, CurvePicker};

	DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, MouseMode& mouse_mode, SelectMode& select_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	~DogEditor();

	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up();

	void render_pairs() const {render_paired_constraints();render_selected_pairs();}
	bool has_constraints() {return (b.rows() + edgePoints.size()) > 0;}

	void render_positional_constraints() const;
	void render_paired_constraints() const;
	void render_selected_pairs() const;
	void clearHandles();
private:
	Dog& dog;
	MouseMode& mouse_mode; SelectMode& select_mode;
	bool& has_new_constraints;
	Eigen::VectorXi& b; Eigen::VectorXd& bc; std::vector<std::pair<int,int>>& paired_vertices;
	std::vector<EdgePoint>& edgePoints; Eigen::MatrixXd& edgeCoords;

	
	
	// TODO: Use V_ren and F_ren with the editor. Convert inner points to the correct V index and edge points to edge point constraints


	void applySelection();
	void onNewHandleID();
	void compute_handle_centroids();
	void get_new_handle_locations();
	Eigen::Vector3f computeTranslation (int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
	Eigen::Vector4f computeRotation(int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);

	igl::opengl::glfw::Viewer& viewer;
	const Eigen::MatrixXd &V; Eigen::MatrixXd oldV;
	Eigen::MatrixXi F; // The triangular faces of the dog (not dogFrendering)
	Lasso lasso;

	Eigen::Vector3f translation;

	//list of currently selected vertices
	Eigen::VectorXi selected_v;//(0,1);
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

	int pair_vertex_1 = -1; int pair_vertex_2 = -1; bool next_pair_first = true;

	//centroids of handle regions, #H x1
	Eigen::MatrixXd handle_centroids;//(0,3);
	int down_mouse_x = -1, down_mouse_y = -1;

	int spline_pt_number = 4;

	bool action_started = false;
};