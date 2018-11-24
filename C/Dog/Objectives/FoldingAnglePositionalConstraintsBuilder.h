#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class FoldingAnglePositionalConstraintsBuilder {
public:
	FoldingAnglePositionalConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, const double& alpha);

	void get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc);
private:
	void setRotAxis(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
				const Eigen::RowVector3d& rotCenter);

	int const_n;
	const double& alpha; // starts at zero and goes up
	// Constrained edge indices
	Edge edge1, edge2;
	double edge_t_coordinate;
	Eigen::RowVector3d edge1_p1,edge1_p2, edge2_p1, edge2_p2; // The initial points
	Eigen::RowVector3d center, axis; // Rotation axis of edge1 and center of rotation
	Eigen::VectorXi b;
};
