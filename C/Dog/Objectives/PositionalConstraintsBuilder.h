#pragma once

#include "../Dog.h"

// This class does positional constraints with time-varying positions
//   linear interpolation from initial position to terminal position
class PositionalConstraintsBuilder {
public:
	PositionalConstraintsBuilder(const double& timestep);

	void add_constraint(int row, const Eigen::RowVector3d& src, const Eigen::RowVector3d& dst);
	void get_constraint_positions(Eigen::VectorXd& bc);
private:
	// between 0 and 1
	const double& timestep;
	std::vector<int> row_i;
	std::vector<Eigen::RowVector3d> src_pos;
	std::vector<Eigen::RowVector3d> dst_pos;
};
