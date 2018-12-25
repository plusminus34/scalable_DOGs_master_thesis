#include "Dog/DogSolver.h"
#include "Gui/Editor.h"

class DeformationController {
public:
	DeformationController(): geoConstraintsBuilder(NULL), dogSolver(NULL) {}
	void single_optimization();
	void init_from_new_dog(igl::opengl::glfw::Viewer& viewer, Dog& dog, const QuadTopology& quadTop);

	enum DeformationType {
		DIHEDRAL_FOLDING = 0,
		CURVE_DEFORMATION = 1
	};

	DeformationController::DeformationType deformationType = CURVE_DEFORMATION;

	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up() {return editor->callback_mouse_up();}

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};
	void update_positional_constraints(bool update_solver = true);
	void render_positional_constraints() const {return editor->render_positional_constraints();}
	void reset_constraints() {editor->clearHandles(); b.resize(0);bc.resize(0); reset_dog_solver();}
	bool has_constraints() {return (b.rows() + edgePoints.size()) > 0;}

	double folding_angle = 0;
	double curve_timestep = 0;

	DogSolver::Params p;
	double constraints_deviation;
	double objective;
	Editor::MouseMode mouse_mode = Editor::NONE;
	Editor::SelectMode select_mode = Editor::VertexPicker;

private:
	void reset_dog_solver();

	igl::opengl::glfw::Viewer* viewer;
	// This needs to be reset when the DOG change, or when the soft positional constraints indices change
	//	Since this amounts to a different objective/hessian sparsity pattern
	// This doesn't change when the values of the soft constraints change
	CurveInterpolationConstraintsBuilder* geoConstraintsBuilder;
	DogSolver* dogSolver;
	Editor* editor;
	//FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;
	//CurveInterpolationConstraintsBuilder curveConstraintsBuilder;

	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};