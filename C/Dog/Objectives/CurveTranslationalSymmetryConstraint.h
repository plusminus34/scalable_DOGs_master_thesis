#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

// Todo:Assume the translation variables are at the end of the mesh for now
class CurveTranslationalSymmetryConstraint : public Constraints {
public:
	CurveTranslationalSymmetryConstraint(std::vector<int>& src_points, std::vector<int>& target_points)
			: src_points(src_points), target_points(target_points) {
		const_n = 3*src_points.size();
		IJV.resize(9*src_points.size());
	};

	virtual CurveTranslationalSymmetryConstraint* clone() const {return new CurveTranslationalSymmetryConstraint(*this);}

	bool is_constraint() {return (src_points.size() > 0);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(const_n); 
		int vnum = (x.rows()-3)/3; // Todo:Assume the translation variables are at the end of the mesh for now
		int tx_i = x.rows()-3,ty_i = tx_i+1, tz_i = ty_i+1;
		double t_x = x(tx_i), t_y = x(ty_i), t_z = x(tz_i); int val_cnt = 0;
		for (int i = 0; i < src_points.size(); i++) {
			vals(val_cnt++) = x(src_points[i])+t_x-x(target_points[i]);
			vals(val_cnt++) = x(src_points[i]+vnum)+t_y-x(target_points[i]+vnum);
			vals(val_cnt++) = x(src_points[i]+2*vnum)+t_z-x(target_points[i]+2*vnum);
		}
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;
		int vnum = (x.rows()-3)/3; // Todo:Assume the translation variables are at the end of the mesh for now
		int tx_i = x.rows()-3,ty_i = tx_i+1, tz_i = ty_i+1;
		double t_x = x(tx_i), t_y = x(ty_i), t_z = x(tz_i);
		for (int i = 0; i < src_points.size(); i++) {
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i], 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i], -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, tx_i, 1);

			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+vnum, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i]+vnum, -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ty_i, 1);

			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+2*vnum, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i]+2*vnum, -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, tz_i, 1);

			const_cnt++;
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};
	std::vector<int> src_points; std::vector<int> target_points;
};