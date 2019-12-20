#include "PositionalConstraintsBuilder.h"

using namespace std;

PositionalConstraintsBuilder::PositionalConstraintsBuilder(const Dog& dog, const double& timestep) :
	timestep(timestep), v_num(dog.get_v_num()) {}

void PositionalConstraintsBuilder::add_constraint(int row, const Eigen::RowVector3d& src, const Eigen::RowVector3d& dst){
	row_i.push_back(row);
	src_pos.push_back(src);
	dst_pos.push_back(dst);
}

void PositionalConstraintsBuilder::get_constraint_positions(Eigen::VectorXd& bc){
	for(int i=0; i<row_i.size(); ++i){
		bc(row_i[i]) = (1.0-timestep)*src_pos[i][0] + timestep*dst_pos[i][0];
		bc(row_i[i] + v_num) = (1.0-timestep)*src_pos[i][1] + timestep*dst_pos[i][1];
		bc(row_i[i] +2*v_num) = (1.0-timestep)*src_pos[i][2] + timestep*dst_pos[i][2];
	}
}
