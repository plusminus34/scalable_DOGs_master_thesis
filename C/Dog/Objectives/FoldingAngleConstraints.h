#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class FoldingAngleConstraints: public Constraints {
public:
	FoldingAngleConstraints(const QuadTopology& quadTop, Eigen::MatrixXd& V, Edge edge1, Edge edge2, std::pair<double,double> edge_coordinates);
	Eigen::VectorXd Vals(const Eigen::VectorXd& x);
	Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x);

	void set_angle(double angle) {alpha = angle;}
private:
	const QuadTopology& quadTop;

	double alpha; // starts at zero and goes up
	// Constrained edge indices
	Edge edge1, edge2;
	std::pair<double,double> edge_coordinates;
	Eigen::RowVector3d edge1_p1,edge1_p2, edge2_p1, edge2_p2; // The initial points
	Eigen::RowVector3d axis, center; // Rotation axis of edge1 and center of rotation
	Eigen::VectorXi b;
};
