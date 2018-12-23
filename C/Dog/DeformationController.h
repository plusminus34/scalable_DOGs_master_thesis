#include "DogSolver.h"

class DeformationController {
public:
	DeformationController(): geoConstraintsBuilder(NULL), dogSolver(NULL) {}
	void single_optimization();
	void init_from_new_dog(Dog& dog, const QuadTopology& quadTop);

	enum DeformationType {
		DIHEDRAL_FOLDING = 0,
		CURVE_DEFORMATION = 1
	};

	DeformationController::DeformationType deformationType = CURVE_DEFORMATION;

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};
	void update_positional_constraints();

	double folding_angle = 0;
	double curve_timestep = 0;

	DogSolver::Params p;
	double constraints_deviation;
	double objective;

private:
	// This needs to be reset when the DOG change, or when the soft positional constraints indices change
	//	Since this amounts to a different objective/hessian sparsity pattern
	// This doesn't change when the values of the soft constraints change
	CurveInterpolationConstraintsBuilder* geoConstraintsBuilder;
	DogSolver* dogSolver;
	//FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;
	//CurveInterpolationConstraintsBuilder curveConstraintsBuilder;

	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};