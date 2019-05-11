#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

// Todo:Assume the translation variables are at the end of the mesh for now
class CurveAffineSymmetryConstraint : public Constraints {
public:
	CurveAffineSymmetryConstraint(std::vector<int>& src_points, std::vector<int>& target_points)
			: src_points(src_points), target_points(target_points) {
		const_n = 3*src_points.size();
		IJV.resize(24*src_points.size());
	};

	virtual CurveAffineSymmetryConstraint* clone() const {return new CurveAffineSymmetryConstraint(*this);}

	bool is_constraint() {return (src_points.size() > 0);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(const_n); 
		int vnum = (x.rows()-12)/3; // Todo:Assume the affine variables are the last 12 for now, and that the translations are the first 3
		int tx_i = x.rows()-12,ty_i = tx_i+1, tz_i = ty_i+1;
		double t_x = x(tx_i), t_y = x(ty_i), t_z = x(tz_i); int val_cnt = 0;
		int r11_i = tz_i+1, r12_i = tz_i+2, r13_i = tz_i+3;
		int r21_i = tz_i+4, r22_i = tz_i+5, r23_i = tz_i+6;
		int r31_i = tz_i+7, r32_i = tz_i+8, r33_i = tz_i+9;

		double r11(x(r11_i)),r12(x(r12_i)),r13(x(r13_i));
		double r21(x(r21_i)),r22(x(r22_i)),r23(x(r23_i));
		double r31(x(r31_i)),r32(x(r32_i)),r33(x(r33_i));
		
		for (int i = 0; i < src_points.size(); i++) {
			double p1_x(x(src_points[i])), p1_y(x(src_points[i]+vnum)),p1_z(x(src_points[i]+2*vnum));
			double p2_x(x(target_points[i])), p2_y(x(target_points[i]+vnum)),p2_z(x(target_points[i]+2*vnum));
			vals(val_cnt++) = -p2_x+t_x+p1_x*r11+p1_y*r21+p1_z*r31;
			vals(val_cnt++) = -p2_y+t_y+p1_x*r12+p1_y*r22+p1_z*r32;
			vals(val_cnt++) = -p2_z+t_z+p1_x*r13+p1_y*r23+p1_z*r33;
		}
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;
		int vnum = (x.rows()-12)/3; // Todo:Assume the affine variables are the last 12 for now, and that the translations are the first 3
		int tx_i = x.rows()-12,ty_i = tx_i+1, tz_i = ty_i+1;
		double t_x = x(tx_i), t_y = x(ty_i), t_z = x(tz_i); int val_cnt = 0;
		int r11_i = tz_i+1, r12_i = tz_i+2, r13_i = tz_i+3;
		int r21_i = tz_i+4, r22_i = tz_i+5, r23_i = tz_i+6;
		int r31_i = tz_i+7, r32_i = tz_i+8, r33_i = tz_i+9;

		double r11(x(r11_i)),r12(x(r12_i)),r13(x(r13_i));
		double r21(x(r21_i)),r22(x(r22_i)),r23(x(r23_i));
		double r31(x(r31_i)),r32(x(r32_i)),r33(x(r33_i));
		for (int i = 0; i < src_points.size(); i++) {
			double p1_x(x(src_points[i])), p1_y(x(src_points[i]+vnum)),p1_z(x(src_points[i]+2*vnum));
			double p2_x(x(target_points[i])), p2_y(x(target_points[i]+vnum)),p2_z(x(target_points[i]+2*vnum));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i], r11);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+vnum, r21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+2*vnum, r31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i], -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, tx_i, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r11_i, p1_x);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r21_i, p1_y);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r31_i, p1_z);

			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i], r12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+vnum, r22);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+2*vnum, r32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i]+vnum, -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ty_i, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r12_i, p1_x);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r22_i, p1_y);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r32_i, p1_z);

			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i], r13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+vnum, r23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, src_points[i]+2*vnum, r33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, target_points[i]+2*vnum, -1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, tz_i, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r13_i, p1_x);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r23_i, p1_y);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, r33_i, p1_z);

			const_cnt++;
		}
	}

	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};
	std::vector<int> src_points; std::vector<int> target_points;
};