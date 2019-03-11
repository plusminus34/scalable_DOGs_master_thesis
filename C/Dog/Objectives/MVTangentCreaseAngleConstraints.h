#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "../../Optimization/Constraints.h"

struct MVTangentCreaseFold {
	int v1,v2;
	EdgePoint ep,ep_b,ep_f;
};

class MVTangentCreaseAngleConstraints : public Constraints {
public:
	MVTangentCreaseAngleConstraints(const Eigen::VectorXd& x, const std::vector<MVTangentCreaseFold>& fold_params, const std::vector<double> cos_angles) :
			fold_params(fold_params), cos_angles(cos_angles) {
		// Every edge pairs consist 1 angle constraint
		const_n = fold_params.size();
		IJV.resize(6*const_n);

		int vnum = x.rows()/3;
		for (int i = 0; i < fold_params.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			EdgePoint ep(fold_params[i].ep), ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double ep_0_t(ep.t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));


			Eigen::RowVector3d v1(v1_x,v1_y,v1_z); Eigen::RowVector3d v2(v2_x,v2_y,v2_z);
			Eigen::RowVector3d ep_0_v = v1*ep_0_t+(1-ep_0_t)*v2;
			Eigen::RowVector3d ep_b_v1(ep_b_v1_x,ep_b_v1_y,ep_b_v1_z); Eigen::RowVector3d ep_b_v2(ep_b_v2_x,ep_b_v2_y,ep_b_v2_z);
			Eigen::RowVector3d ep_b_v = ep_b_v1*ep_b_t+(1-ep_0_t)*ep_b_v2;
			Eigen::RowVector3d ep_f_v1(ep_f_v1_x,ep_f_v1_y,ep_f_v1_z); Eigen::RowVector3d ep_f_v2(ep_f_v2_x,ep_f_v2_y,ep_f_v2_z);
			Eigen::RowVector3d ep_f_v = ep_f_v1*ep_f_t+(1-ep_0_t)*ep_f_v2;


			double l1((ep_0_v-ep_b_v).norm()),l2((ep_f_v-ep_0_v).norm());
			curve_l1.push_back(l1); curve_l2.push_back(l2);
			double fold_e_l = (v1-v2).norm();
			e_lens.push_back(fold_e_l);

			auto curve_T = (l2*(ep_0_v-ep_b_v)+l1*(ep_f_v-ep_0_v))/(l1*l2);
			double dot_T_grid_t = curve_T.dot(v1-v2);
			double curve_T_grid_t_angle = acos(dot_T_grid_t/fold_e_l); // The other T is already normalized
			fold_e_crease_angles.push_back(curve_T_grid_t_angle);
		};
	};

	void set_angles(const std::vector<double> cos_angles_i) { cos_angles = cos_angles_i;}
	virtual MVTangentCreaseAngleConstraints* clone() const {return new MVTangentCreaseAngleConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		int vnum = x.rows()/3;
		Eigen::VectorXd vals(fold_params.size()); 
		for (int i = 0; i < edge_pairs.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			EdgePoint ep(fold_params[i].ep), ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double ep_0_t(ep.t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double cos_angle = cos_angles[i];
			double fold_e_l = e_lens[i];
			double l1 = curve_l1[i]; double l2 = curve_l2[i];
			double fold_e_crease_angle = fold_e_crease_angles[i];

			double t2 = 1.0/fold_e_l;
			double t3 = ep_0_t-1.0;
			double t4 = t3*v2_x;
			double t5 = 1.0/l1;
			double t6 = 1.0/l2;
			double t7 = ep_b_t-1.0;
			double t8 = ep_f_t-1.0;
			double t9 = t3*v2_y;
			double t10 = sin(fold_e_crease_angle);
			double t11 = 1.0/t10;
			double t12 = v1_y-v2_y;
			double t13 = v1_z-v2_z;
			double t14 = ep_b_t*ep_b_v1_x;
			double t17 = ep_0_t*v1_x;
			double t15 = t4+t14-t17-ep_b_v2_x*t7;
			double t16 = ep_f_t*ep_f_v1_x;
			double t18 = l2*t15;
			double t19 = v1_x-v2_x;
			double t20 = t3*v2_z;
			double t21 = ep_b_t*ep_b_v1_y;
			double t24 = ep_0_t*v1_y;
			double t22 = t9+t21-t24-ep_b_v2_y*t7;
			double t23 = ep_f_t*ep_f_v1_y;
			double t25 = l2*t22;
			double t26 = ep_b_t*ep_b_v1_z;
			double t29 = ep_0_t*v1_z;
			double t27 = t20+t26-t29-ep_b_v2_z*t7;
			double t28 = ep_f_t*ep_f_v1_z;
			vals(i) = -cos_angle+t2*t11*t19*(t2*t5*t6*t13*(t25-l1*(t9+t23-t24-ep_f_v2_y*t8))+t2*t5*t6*t12*(l1*(t20+t28-t29-ep_f_v2_z*t8)-l2*t27))+t2*t11*t13*(t2*t5*t6*t12*(t18-l1*(t4+t16-ep_f_v2_x*t8-ep_0_t*v1_x))-t2*t5*t6*t19*(t25-l1*(t9+t23-ep_f_v2_y*t8-ep_0_t*v1_y)))-t2*t11*t12*(t2*t5*t6*t13*(t18-l1*(t4+t16-t17-ep_f_v2_x*t8))-t2*t5*t6*t19*(l2*t27-l1*(t20+t28-ep_f_v2_z*t8-ep_0_t*v1_z)));

			//std::cout << "dot prod = " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z)) << std::endl;
			//std::cout << "cos_angle = " << cos_angle << std::endl;
			//std::cout << "normalized it's: " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z))/e_len << std::endl;
			//exit(1);
		};
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;

		int vnum = x.rows()/3;
		Eigen::VectorXd vals(edge_pairs.size()); 
		for (int i = 0; i < edge_pairs.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			EdgePoint ep(fold_params[i].ep), ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double ep_0_t(ep.t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double cos_angle = cos_angles[i];
			double fold_e_l = e_lens[i];
			double l1 = curve_l1[i]; double l2 = curve_l2[i];
			double fold_e_crease_angle = fold_e_crease_angles[i];

			double t2 = 1.0/fold_e_l;
			double t3 = ep_0_t-1.0;
			double t4 = t3*v2_y;
			double t5 = 1.0/l1;
			double t6 = 1.0/l2;
			double t7 = ep_b_t-1.0;
			double t8 = ep_f_t-1.0;
			double t9 = t3*v2_z;
			double t10 = sin(fold_e_crease_angle);
			double t11 = 1.0/t10;
			double t12 = ep_b_t*ep_b_v1_y;
			double t15 = ep_0_t*v1_y;
			double t44 = ep_b_v2_y*t7;
			double t13 = t4+t12-t15-t44;
			double t14 = ep_f_t*ep_f_v1_y;
			double t16 = l2*t13;
			double t17 = v1_y-v2_y;
			double t18 = v1_z-v2_z;
			double t19 = ep_b_t*ep_b_v1_z;
			double t22 = ep_0_t*v1_z;
			double t27 = ep_b_v2_z*t7;
			double t20 = t9+t19-t22-t27;
			double t21 = ep_f_t*ep_f_v1_z;
			double t23 = l2*t20;
			double t24 = ep_0_t*l1;
			double t36 = ep_0_t*l2;
			double t25 = t24-t36;
			double t26 = t3*v2_x;
			double t37 = ep_f_v2_z*t8;
			double t28 = t9+t21-t22-t37;
			double t38 = l1*t28;
			double t29 = t23-t38;
			double t30 = ep_b_t*ep_b_v1_x;
			double t33 = ep_0_t*v1_x;
			double t41 = ep_b_v2_x*t7;
			double t31 = t26+t30-t33-t41;
			double t32 = ep_f_t*ep_f_v1_x;
			double t34 = l2*t31;
			double t35 = v1_x-v2_x;
			double t39 = t2*t5*t6*t29;
			double t40 = t39-t2*t5*t6*t18*t25;
			double t47 = ep_f_v2_x*t8;
			double t42 = t26+t32-t33-t47;
			double t48 = l1*t42;
			double t43 = t34-t48;
			double t51 = ep_f_v2_y*t8;
			double t45 = t4+t14-t15-t51;
			double t52 = l1*t45;
			double t46 = t16-t52;
			double t49 = t2*t5*t6*t43;
			double t50 = t49-t2*t5*t6*t25*t35;
			double t53 = t2*t5*t6*t46;
			double t54 = t53-t2*t5*t6*t17*t25;
			double t55 = l1*t3;
			double t57 = l2*t3;
			double t56 = t55-t57;
			double t58 = t39-t2*t5*t6*t18*t56;
			double t59 = t2*t5*t6*t17*t43;
			double t60 = t59-t2*t5*t6*t35*t46;
			double t61 = t2*t11*t60;
			double t62 = t49-t2*t5*t6*t35*t56;
			double t63 = t53-t2*t5*t6*t17*t56;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i, t2*t11*(t2*t5*t6*t18*(t16-l1*(t4+t14-ep_f_v2_y*t8-ep_0_t*v1_y))-t2*t5*t6*t17*(t23-l1*(t9+t21-ep_f_v2_z*t8-ep_0_t*v1_z)))+t2*t11*t17*t40-t2*t11*t18*t54);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+vnum, -t2*t11*(t2*t5*t6*t18*(t34-l1*(t26+t32-ep_f_v2_x*t8-ep_0_t*v1_x))-t2*t5*t6*t29*t35)+t2*t11*t18*t50-t2*t11*t35*t40);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+2*vnum, t61-t2*t11*t17*t50+t2*t11*t35*t54);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i, t2*t11*(t2*t5*t6*t17*t29-t2*t5*t6*t18*t46)-t2*t11*t17*t58+t2*t11*t18*t63);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+vnum, t2*t11*(t2*t5*t6*t18*t43-t2*t5*t6*t29*t35)-t2*t11*t18*t62+t2*t11*t35*t58);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+2*vnum, -t61+t2*t11*t17*t62-t2*t11*t35*t63);

  			const_cnt++;
  		}
  		if (const_cnt != const_n) {
			std::cout << "MVTangentCreaseAngleConstraints: Error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << std::endl;
			exit(1);
		}
		if (ijv_cnt != IJV.size()) {
			std::cout << "MVTangentCreaseAngleConstraints: Error! ijv_cnt = " << ijv_cnt << " but IJV.size() " << IJV.size() << std::endl;
			exit(1);
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// The constraints are not linear but for now we use here gauss-newton..
	};
	const std::vector<MVTangentCreaseFold>& fold_params;
	std::vector<std::pair<Edge,Edge>> edge_pairs;
	std::vector<double> cos_angles;
	std::vector<double> e_lens;
	std::vector<double> curve_l1,curve_l2;
	std::vector<double> fold_e_crease_angles;
};