#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

// TODO: Use positional constraints (update b every time and have positioal constraints)
// Add that to set_angle

class FoldingAngleConstraints: public Constraints {
public:
	FoldingAngleConstraints(const Eigen::MatrixXd& V, Edge edge1, Edge edge2, double edge_t_coordinate);
	Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const;

	void set_angle(double angle) {alpha = angle;}
private:
	static Eigen::RowVector3d rotate_vec(const Eigen::RowVector3d& pt, const Eigen::RowVector3d& center, const Eigen::Vector3d& axis, double angle);

	double alpha; // starts at zero and goes up
	// Constrained edge indices
	Edge edge1, edge2;
	double edge_t_coordinate;
	Eigen::RowVector3d edge1_p1,edge1_p2, edge2_p1, edge2_p2; // The initial points
	Eigen::RowVector3d axis, center; // Rotation axis of edge1 and center of rotation
	Eigen::VectorXi b;
};
