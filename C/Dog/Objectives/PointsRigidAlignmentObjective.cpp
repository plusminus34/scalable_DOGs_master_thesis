#include "PointsRigidAlignmentObjective.h"
#include <igl/procrustes.h>


PointsRigidAlignmentObjective::PointsRigidAlignmentObjective(std::vector<int>& src_points, std::vector<int>& target_points) :
		src_points(src_points), target_points(target_points) { 

	IJV.resize(src_points.size()*174);
}

void PointsRigidAlignmentObjective::update_rigid_motion(const Eigen::VectorXd& x, 
			Eigen::Matrix3d& R, Eigen::Vector3d& T) const {
	int vnum = x.rows()/3;
	Eigen::MatrixXd src(src_points.size(),3); Eigen::MatrixXd target(target_points.size(),3);
	for (int i = 0; i < src_points.size(); i++) {
		src.row(i) << x(src_points[i]), x(src_points[i]+vnum), x(src_points[i]+2*vnum);
		target.row(i) << x(target_points[i]), x(target_points[i]+vnum), x(target_points[i]+2*vnum);
	}
	
	// call procrustes with includeScaling = false, includeReflections = false
	double scale_dummy; igl::procrustes(src, target, false, false, scale_dummy,R,T);
}

double PointsRigidAlignmentObjective::obj(const Eigen::VectorXd& x) const {
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,R,T);
	return 0;
}

Eigen::VectorXd PointsRigidAlignmentObjective::grad(const Eigen::VectorXd& x) const {
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,R,T);
	Eigen::VectorXd grad;
  	grad.resize(x.rows(),1); grad.setZero();
  	int vnum = x.rows()/3;

  	return grad;
}


void PointsRigidAlignmentObjective::updateHessianIJV(const Eigen::VectorXd& x) {
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,R,T);
}
