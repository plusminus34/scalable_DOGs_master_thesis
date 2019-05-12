#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

// Todo:Assume the translation variables are at the end of the mesh for now
class Curve2AffineSymmetriesConstraint : public Constraints {
public:
	Curve2AffineSymmetriesConstraint(std::vector<int>& src_points1, std::vector<int>& target_points1,
									std::vector<int>& src_points2, std::vector<int>& target_points2)
			: src_points1(src_points1), target_points1(target_points1), src_points2(src_points2), target_points2(target_points2) {
		const_n = 3*(src_points1.size()+src_points2.size())+12; // Including affine commuting constraints
		IJV.resize(24*(src_points1.size()+src_points2.size())+120); // Including affine commuting constraints
	};

	virtual Curve2AffineSymmetriesConstraint* clone() const {return new Curve2AffineSymmetriesConstraint(*this);}

	bool is_constraint() {return (src_points1.size() > 0);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(const_n); 
		// Assume the affine variables are the last 24 for now, and that the translations are the first 3
		// And that we have all the first affine variables before the second
		int vnum = (x.rows()-24)/3; int val_cnt = 0;
		set_affine_vals(vals, x, src_points1, target_points1, x.rows()-24, vnum, val_cnt);
		set_affine_vals(vals, x, src_points2, target_points2, x.rows()-12, vnum, val_cnt);

		
		// add affine motion commuting constraints
		{
			int t1_1i = x.rows()-24,t1_2i = t1_1i+1, t1_3i = t1_2i+1;
			double t1_1 = x(t1_1i), t1_2 = x(t1_2i), t1_3 = x(t1_2i);
			int r1_11i = t1_3i+1, r1_12i = t1_3i+2, r1_13i = t1_3i+3;
			int r1_21i = t1_3i+4, r1_22i = t1_3i+5, r1_23i = t1_3i+6;
			int r1_31i = t1_3i+7, r1_32i = t1_3i+8, r1_33i = t1_3i+9;

			double r1_11(x(r1_11i)),r1_12(x(r1_12i)),r1_13(x(r1_13i));
			double r1_21(x(r1_21i)),r1_22(x(r1_22i)),r1_23(x(r1_23i));
			double r1_31(x(r1_31i)),r1_32(x(r1_32i)),r1_33(x(r1_33i));

			int t2_1i = x.rows()-24,t2_2i = t2_1i+1, t2_3i = t2_2i+1;
			double t2_1 = x(t2_1i), t2_2 = x(t2_2i), t2_3 = x(t2_2i);
			int r2_11i = t2_3i+1, r2_12i = t2_3i+2, r2_13i = t2_3i+3;
			int r2_21i = t2_3i+4, r2_22i = t2_3i+5, r2_23i = t2_3i+6;
			int r2_31i = t2_3i+7, r2_32i = t2_3i+8, r2_33i = t2_3i+9;

			double r2_11(x(r2_11i)),r2_12(x(r2_12i)),r2_13(x(r2_13i));
			double r2_21(x(r2_21i)),r2_22(x(r2_22i)),r2_23(x(r2_23i));
			double r2_31(x(r2_31i)),r2_32(x(r2_32i)),r2_33(x(r2_33i));

			double t2 = r1_12*r2_21;
			double t3 = r1_13*r2_31;
			double t4 = r1_23*r2_32;

			vals(val_cnt++) = t2+t3-r1_21*r2_12-r1_31*r2_13;
			vals(val_cnt++) = -r1_11*r2_21+r1_21*r2_11-r1_21*r2_22+r1_22*r2_21+r1_23*r2_31-r1_31*r2_23;
			vals(val_cnt++) = -r1_11*r2_31+r1_31*r2_11-r1_21*r2_32+r1_32*r2_21-r1_31*r2_33+r1_33*r2_31;
			vals(val_cnt++) = r1_11*r2_12-r1_12*r2_11+r1_12*r2_22-r1_22*r2_12+r1_13*r2_32-r1_32*r2_13;
			vals(val_cnt++) = -t2+t4+r1_21*r2_12-r1_32*r2_23;
			vals(val_cnt++) = -r1_12*r2_31+r1_31*r2_12-r1_22*r2_32+r1_32*r2_22-r1_32*r2_33+r1_33*r2_32;
			vals(val_cnt++) = r1_11*r2_13-r1_13*r2_11+r1_12*r2_23-r1_23*r2_12+r1_13*r2_33-r1_33*r2_13;
			vals(val_cnt++) = -r1_13*r2_21+r1_21*r2_13+r1_22*r2_23-r1_23*r2_22+r1_23*r2_33-r1_33*r2_23;
			vals(val_cnt++) = -t3-t4+r1_31*r2_13+r1_32*r2_23;
			vals(val_cnt++) = t1_1-t2_1+r1_11*t2_1+r1_21*t2_2+r1_31*t2_3-r2_11*t1_1-r2_21*t1_2-r2_31*t1_3;
			vals(val_cnt++) = t1_2-t2_2+r1_12*t2_1+r1_22*t2_2+r1_32*t2_3-r2_12*t1_1-r2_22*t1_2-r2_32*t1_3;
			vals(val_cnt++) = t1_3-t2_3+r1_13*t2_1+r1_23*t2_2+r1_33*t2_3-r2_13*t1_1-r2_23*t1_2-r2_33*t1_3;
		}
		std::cout << "vals.head(12) = " << vals.head(12) << std::endl;
		
		return vals;
	}

	void set_affine_vals(Eigen::VectorXd& vals, const Eigen::VectorXd& x, 
						const std::vector<int>& src_points, const std::vector<int>& target_points, int affine_vars_offset, 
						int vnum, int& val_cnt) const {
		int tx_i = affine_vars_offset,ty_i = tx_i+1, tz_i = ty_i+1;
		double t_x = x(tx_i), t_y = x(ty_i), t_z = x(tz_i);
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
	}

	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;
		int vnum = (x.rows()-24)/3;
		set_affine_jacobian_ijv(ijv_cnt, const_cnt, x, src_points1, target_points1, x.rows()-24, vnum);
		set_affine_jacobian_ijv(ijv_cnt, const_cnt, x, src_points2, target_points2, x.rows()-12, vnum);

		// add affine motion commuting constraints
		{
			int t1_1i = x.rows()-24,t1_2i = t1_1i+1, t1_3i = t1_2i+1;
			double t1_1 = x(t1_1i), t1_2 = x(t1_2i), t1_3 = x(t1_2i);
			int r1_11i = t1_3i+1, r1_12i = t1_3i+2, r1_13i = t1_3i+3;
			int r1_21i = t1_3i+4, r1_22i = t1_3i+5, r1_23i = t1_3i+6;
			int r1_31i = t1_3i+7, r1_32i = t1_3i+8, r1_33i = t1_3i+9;

			double r1_11(x(r1_11i)),r1_12(x(r1_12i)),r1_13(x(r1_13i));
			double r1_21(x(r1_21i)),r1_22(x(r1_22i)),r1_23(x(r1_23i));
			double r1_31(x(r1_31i)),r1_32(x(r1_32i)),r1_33(x(r1_33i));

			int t2_1i = x.rows()-24,t2_2i = t2_1i+1, t2_3i = t2_2i+1;
			double t2_1 = x(t2_1i), t2_2 = x(t2_2i), t2_3 = x(t2_2i);
			int r2_11i = t2_3i+1, r2_12i = t2_3i+2, r2_13i = t2_3i+3;
			int r2_21i = t2_3i+4, r2_22i = t2_3i+5, r2_23i = t2_3i+6;
			int r2_31i = t2_3i+7, r2_32i = t2_3i+8, r2_33i = t2_3i+9;

			double r2_11(x(r2_11i)),r2_12(x(r2_12i)),r2_13(x(r2_13i));
			double r2_21(x(r2_21i)),r2_22(x(r2_22i)),r2_23(x(r2_23i));
			double r2_31(x(r2_31i)),r2_32(x(r2_32i)),r2_33(x(r2_33i));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, -r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, -r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, -r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, -r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, r1_13);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_11i, -r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, r2_11-r2_22);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_22i, r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, -r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_11i, r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, -r1_11+r1_22);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_22i, -r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, -r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, r1_23);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_11i, -r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, -r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, r2_11-r2_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_33i, r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_11i, r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, -r1_11+r1_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, -r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_33i, -r1_31);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_11i, r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, -r2_11+r2_22);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_22i, -r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, -r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_11i, -r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, r1_11-r1_22);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, -r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_22i, r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, r1_13);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, -r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, -r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, -r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, -r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, r1_23);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, -r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_22i, -r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, r2_22-r2_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_33i, r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_22i, r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, -r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, -r1_22+r1_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_33i, -r1_32);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_11i, r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, -r2_11+r2_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, -r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_33i, -r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_11i, -r1_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, -r1_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, r1_11-r1_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_33i, r1_13);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, -r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_22i, r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, -r2_22+r2_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_33i, -r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, -r1_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_22i, -r1_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, r1_22-r1_33);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_33i, r1_23);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, -r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, -r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, -r1_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, -r1_23);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_1i, -r2_11+1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_2i, -r2_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_3i, -r2_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_11i, t2_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_21i, t2_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_31i, t2_3);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_1i, r1_11-1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_2i, r1_21);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_3i, r1_31);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_11i, -t1_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_21i, -t1_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_31i, -t1_3);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_1i, -r2_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_2i, -r2_22+1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_3i, -r2_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_12i, t2_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_22i, t2_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_32i, t2_3);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_1i, r1_12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_2i, r1_22-1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_3i, r1_32);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_12i, -t1_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_22i, -t1_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_32i, -t1_3);
			const_cnt++;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_1i, -r2_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_2i, -r2_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t1_3i, -r2_33+1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_13i, t2_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_23i, t2_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r1_33i, t2_3);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_1i, r1_13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_2i, r1_23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,t2_3i, r1_33-1.0);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_13i, -t1_1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_23i, -t1_2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,r2_33i, -t1_3);
			const_cnt++;

		}

	}

	void set_affine_jacobian_ijv(int& ijv_cnt, int& const_cnt, const Eigen::VectorXd& x, 
						const std::vector<int>& src_points, const std::vector<int>& target_points, int affine_vars_offset, int vnum) {

		int tx_i = affine_vars_offset,ty_i = tx_i+1, tz_i = ty_i+1;
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
	std::vector<int> src_points1; std::vector<int> target_points1;
	std::vector<int> src_points2; std::vector<int> target_points2;
};