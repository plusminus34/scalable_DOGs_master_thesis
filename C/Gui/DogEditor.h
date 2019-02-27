#pragma once

#include "../Dog/DogSolver.h"
#include "../Folding/CurvedFoldingBiasObjective.h"
#include "Editor.h"

class DogEditor {
public:
	DogEditor(): dogSolver(NULL) {}
	~DogEditor();
	void single_optimization();
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}
	void init_from_new_dog(Dog& dog);

	enum DeformationType {
		DIHEDRAL_FOLDING = 0,
		CURVE_DEFORMATION = 1
	};

	DogEditor::DeformationType deformationType = CURVE_DEFORMATION;

	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up() {return editor->callback_mouse_up();}

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};
	void render_positional_constraints() const {return editor->render_positional_constraints();}
	void render_pairs() const {editor->render_paired_constraints();editor->render_selected_pairs();}
	void reset_constraints() {b.resize(0);bc.resize(0); paired_vertices.clear(); edgePoints.clear(); edgeCoords.resize(0,3); editor->clearHandles(); reset_dog_solver();}
	bool has_constraints() {return (b.rows() + edgePoints.size()) > 0;}

	void add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc);
	void add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords);
	void add_edge_point_constraint(const EdgePoint& new_edgePoints, const Eigen::RowVector3d& new_edgeCoords);
	void add_pair_vertices_constraints(const std::vector<std::pair<int,int>>& new_pair_vertices);
	void add_pair_vertices_constraint(int v1, int v2);

	void update_edge_coords(Eigen::MatrixXd& edgeCoords_i) {edgeCoords = edgeCoords_i; dogSolver->update_edge_coords(edgeCoords);}
	void update_point_coords(Eigen::VectorXd& bc_i) {bc = bc_i; dogSolver->update_point_coords(bc);}

	double curve_timestep = 0;

	DogSolver::Params p;
	double constraints_deviation;
	double objective;
	Editor::MouseMode mouse_mode = Editor::NONE;
	Editor::SelectMode select_mode = Editor::VertexPicker;

private:
	void reset_dog_solver();

	igl::opengl::glfw::Viewer* viewer;
	Eigen::MatrixXi FTriangular; // The triangular faces of the dog (not dogFrendering)
	// TODO: Use V_ren and F_ren with the editor. Convert inner points to the correct V index and edge points to edge point constraints
	// This needs to be reset when the DOG change, or when the soft positional constraints indices change
	//	Since this amounts to a different objective/hessian sparsity pattern
	// This doesn't change when the values of the soft constraints change
	DogSolver* dogSolver;
	Editor* editor;

	Eigen::VectorXd init_x0;
	//FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;
	//CurveInterpolationConstraintsBuilder curveConstraintsBuilder;

	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	// Point pair constraints
	std::vector<std::pair<int,int>> paired_vertices;
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};