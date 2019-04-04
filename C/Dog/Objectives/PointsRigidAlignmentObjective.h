#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class PointsRigidAlignmentObjective : public Objective {
  
public:
	PointsRigidAlignmentObjective(std::vector<int>& src_points, std::vector<int>& target_points);
	virtual PointsRigidAlignmentObjective* clone() const {return new PointsRigidAlignmentObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

	static void update_rigid_motion(const Eigen::VectorXd& x, const std::vector<int>& src_points, const std::vector<int>& target_points,
				Eigen::Matrix3d& R, Eigen::Vector3d& T);
private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);
	
	std::vector<int> src_points; std::vector<int> target_points;
};