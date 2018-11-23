#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class FoldingAnglePositionalConstraintsBuilder {
public:
	FoldingAnglePositionalConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS);

	void get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc);

	void set_angle(double angle) {alpha = angle;}
private:
	static Eigen::RowVector3d rotate_vec(const Eigen::RowVector3d& pt, const Eigen::RowVector3d& center, const Eigen::Vector3d& axis, double angle);
	void setRotAxis(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
				const Eigen::RowVector3d& rotCenter);

	int const_n;
	double alpha; // starts at zero and goes up
	// Constrained edge indices
	Edge edge1, edge2;
	double edge_t_coordinate;
	Eigen::RowVector3d edge1_p1,edge1_p2, edge2_p1, edge2_p2; // The initial points
	Eigen::RowVector3d center, axis; // Rotation axis of edge1 and center of rotation
	Eigen::VectorXi b;
};
