#pragma once

#include "../Dog/DogSolver.h"
#include "Editor.h"

class DogEditor {
public:
	DogEditor(Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	~DogEditor();

	void init_from_new_dog(Dog& dog);
	bool has_new_constraints() {return editor->new_constraints;}
	void signal_handled_new_constraints() {editor->new_constraints = false;}
	bool has_new_point_constraints() {return new_point_constraints;}

	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}


	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up() {return editor->callback_mouse_up();}

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};
	void render_positional_constraints() const {return editor->render_positional_constraints();}
	void render_pairs() const {editor->render_paired_constraints();editor->render_selected_pairs();}
	bool has_constraints() {return (b.rows() + edgePoints.size()) > 0;}

	void add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords);
	/*
	void add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc);
	
	void add_edge_point_constraint(const EdgePoint& new_edgePoints, const Eigen::RowVector3d& new_edgeCoords);
	void add_pair_vertices_constraints(const std::vector<std::pair<int,int>>& new_pair_vertices);
	void add_pair_vertices_constraint(int v1, int v2);
	*/

	Editor::MouseMode mouse_mode = Editor::NONE;
	Editor::SelectMode select_mode = Editor::VertexPicker;
	Editor* editor;


private:
	Eigen::VectorXi& b; Eigen::VectorXd& bc; std::vector<std::pair<int,int>>& paired_vertices;
	std::vector<EdgePoint>& edgePoints; Eigen::MatrixXd& edgeCoords;

	igl::opengl::glfw::Viewer* viewer;
	Eigen::MatrixXi FTriangular; // The triangular faces of the dog (not dogFrendering)
	// TODO: Use V_ren and F_ren with the editor. Convert inner points to the correct V index and edge points to edge point constraints
	
	bool new_point_constraints;

};