#pragma once

#include "../Dog/DogSolver.h"
#include "Editor.h"

class DogEditor {
public:
	DogEditor(bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
				std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	~DogEditor();

	void init_from_new_dog(Dog& dog);

	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}


	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up() {return editor->callback_mouse_up();}

	void render_positional_constraints() const {return editor->render_positional_constraints();}
	void render_pairs() const {editor->render_paired_constraints();editor->render_selected_pairs();}
	bool has_constraints() {return (b.rows() + edgePoints.size()) > 0;}

	Editor::MouseMode mouse_mode = Editor::NONE;
	Editor::SelectMode select_mode = Editor::VertexPicker;
	Editor* editor;
private:
	bool& has_new_constraints;
	Eigen::VectorXi& b; Eigen::VectorXd& bc; std::vector<std::pair<int,int>>& paired_vertices;
	std::vector<EdgePoint>& edgePoints; Eigen::MatrixXd& edgeCoords;

	igl::opengl::glfw::Viewer* viewer;
	Eigen::MatrixXi FTriangular; // The triangular faces of the dog (not dogFrendering)
	// TODO: Use V_ren and F_ren with the editor. Convert inner points to the correct V index and edge points to edge point constraints
};