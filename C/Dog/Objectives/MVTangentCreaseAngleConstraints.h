#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "../../Optimization/Constraints.h"

struct MVTangentCreaseFold {
	MVTangentCreaseFold(int v1, int v2, int w1, int w2, double edge_t, EdgePoint ep_b, EdgePoint ep_f) : 
			v1(v1),v2(v2),w1(w1),w2(w2),edge_t(edge_t),ep_b(ep_b),ep_f(ep_f) {}
	int v1,v2;
	int w1,w2;
	double edge_t;
	EdgePoint ep_b,ep_f;
};

class MVTangentCreaseAngleConstraints : public Constraints {
public:
	MVTangentCreaseAngleConstraints(const Eigen::VectorXd& x, const std::vector<MVTangentCreaseFold>& fold_params, const std::vector<double> cos_angles) :
			fold_params(fold_params), cos_angles(cos_angles) {
		// Every edge pairs consist 1 angle constraint
		const_n = fold_params.size();
		IJV.resize(24*const_n);

		int vnum = x.rows()/3;
		for (int i = 0; i < fold_params.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			int w1_i(fold_params[i].w1),w2_i(fold_params[i].w2);
			EdgePoint ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));
			const double ep_0_t(fold_params[i].edge_t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));


			Eigen::RowVector3d v1(v1_x,v1_y,v1_z); Eigen::RowVector3d v2(v2_x,v2_y,v2_z);
			Eigen::RowVector3d w1(w1_x,w1_y,w1_z); Eigen::RowVector3d w2(w2_x,w2_y,w2_z);
			Eigen::RowVector3d ep_0_v = v1*ep_0_t+(1-ep_0_t)*v2;
			Eigen::RowVector3d ep_b_v1(ep_b_v1_x,ep_b_v1_y,ep_b_v1_z); Eigen::RowVector3d ep_b_v2(ep_b_v2_x,ep_b_v2_y,ep_b_v2_z);
			Eigen::RowVector3d ep_b_v = ep_b_v1*ep_b_t+(1-ep_0_t)*ep_b_v2;
			Eigen::RowVector3d ep_f_v1(ep_f_v1_x,ep_f_v1_y,ep_f_v1_z); Eigen::RowVector3d ep_f_v2(ep_f_v2_x,ep_f_v2_y,ep_f_v2_z);
			Eigen::RowVector3d ep_f_v = ep_f_v1*ep_f_t+(1-ep_0_t)*ep_f_v2;

			double l1((ep_0_v-ep_b_v).norm()),l2((ep_f_v-ep_0_v).norm());
			curve_l1.push_back(l1); curve_l2.push_back(l2);
			double fold_e1_l = (v1-v2).norm(), fold_e2_l = (w1-w2).norm();
			folds_e1_l.push_back(fold_e1_l); folds_e2_l.push_back(fold_e2_l);

			Eigen::RowVector3d curve_T = (l2*(ep_0_v-ep_b_v)+l1*(ep_f_v-ep_0_v)).normalized();
			double dot_T_grid_t = curve_T.dot(v1-v2);
			//std::cout << "dot_T_grid_t = " << dot_T_grid_t << std::endl;
			double curve_T_grid_t_angle = acos(dot_T_grid_t/fold_e1_l); // The other T is already normalized
			//std::cout << "curve_T_grid_t_angle = " << curve_T_grid_t_angle << std::endl;
			fold_e_crease_angles.push_back(curve_T_grid_t_angle);

			/*
			std::cout << "ep_b_v1 = " << ep_b_v1 << std::endl;
			std::cout << "ep_b_v2 = " << ep_b_v2 << std::endl;

			std::cout << "v1 = " << v1 << std::endl;
			std::cout << "v2 = " << v2 << std::endl;

			std::cout << "w1 = " << v1 << std::endl;
			std::cout << "w2 = " << v2 << std::endl;

			std::cout << "ep_f_v1 = " << ep_f_v1 << std::endl;
			std::cout << "ep_f_v2 = " << ep_f_v2 << std::endl;

			std::cout << "ep_b_t = " << ep_b_t << std::endl;
			std::cout << "ep_f_t = " << ep_f_t << std::endl;

			std::cout << "fold_e_crease_angle = " << curve_T_grid_t_angle << std::endl;
			//std::cout << "cos_angle = " << cos_angle << std::endl;
			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << std::endl;
			std::cout << "fold_e_1 = " << fold_e1_l << std::endl;
			std::cout << "fold_e_2 = " << fold_e2_l << std::endl;
			*/
			//int wait; std::cin >> wait;
		};
	};

	void set_angles(const std::vector<double> cos_angles_i) { cos_angles = cos_angles_i;}
	virtual MVTangentCreaseAngleConstraints* clone() const {return new MVTangentCreaseAngleConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		int vnum = x.rows()/3;
		Eigen::VectorXd vals(fold_params.size()); 
		for (int i = 0; i < fold_params.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			int w1_i(fold_params[i].w1),w2_i(fold_params[i].w2);
			EdgePoint ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));
			const double ep_0_t(fold_params[i].edge_t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double cos_angle = cos_angles[i];
			double fold_e_1 = folds_e1_l[i], fold_e_2 = folds_e2_l[i];
			double l1 = curve_l1[i]; double l2 = curve_l2[i];
			double fold_e_crease_angle = fold_e_crease_angles[i];

			/*
			std::cout << "cos_angle = " << cos_angle << std::endl;
			std::cout << "fold_e_1 = " << fold_e_1 << std::endl;
			std::cout << "fold_e_2 = " << fold_e_2 << std::endl;
			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << std::endl;
			std::cout << "fold_e_crease_angle = " << fold_e_crease_angle << std::endl;
			*/

			double t2 = ep_0_t-1.0;
			double t3 = t2*v2_x;
			double t4 = 1.0/fold_e_2;
			double t5 = 1.0/l1;
			double t6 = 1.0/l2;
			double t7 = ep_b_t-1.0;
			double t8 = ep_f_t-1.0;
			double t9 = t2*v2_y;
			double t10 = 1.0/fold_e_1;
			double t11 = sin(fold_e_crease_angle);
			double t12 = 1.0/t11;
			double t13 = ep_b_t*ep_b_v1_x;
			double t16 = ep_0_t*v1_x;
			double t14 = t3+t13-t16-ep_b_v2_x*t7;
			double t15 = ep_f_t*ep_f_v1_x;
			double t17 = l2*t14;
			double t18 = w1_x-w2_x;
			double t19 = t2*v2_z;
			double t20 = w1_z-w2_z;
			double t21 = ep_b_t*ep_b_v1_y;
			double t24 = ep_0_t*v1_y;
			double t22 = t9+t21-t24-ep_b_v2_y*t7;
			double t23 = ep_f_t*ep_f_v1_y;
			double t25 = l2*t22;
			double t26 = w1_y-w2_y;
			double t27 = ep_b_t*ep_b_v1_z;
			double t30 = ep_0_t*v1_z;
			double t28 = t19+t27-t30-ep_b_v2_z*t7;
			double t29 = ep_f_t*ep_f_v1_z;
			vals(i) = -cos_angle+t10*t12*(t4*t5*t6*t20*(t25-l1*(t9+t23-t24-ep_f_v2_y*t8))+t4*t5*t6*t26*(l1*(t19+t29-t30-ep_f_v2_z*t8)-l2*t28))*(v1_x-v2_x)+t10*t12*(t4*t5*t6*t26*(t17-l1*(t3+t15-ep_f_v2_x*t8-ep_0_t*v1_x))-t4*t5*t6*t18*(t25-l1*(t9+t23-ep_f_v2_y*t8-ep_0_t*v1_y)))*(v1_z-v2_z)-t10*t12*(v1_y-v2_y)*(t4*t5*t6*t20*(t17-l1*(t3+t15-t16-ep_f_v2_x*t8))-t4*t5*t6*t18*(l2*t28-l1*(t19+t29-ep_f_v2_z*t8-ep_0_t*v1_z)));

			//std::cout << "dot prod = " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z)) << std::endl;
			std::cout << "dest cos angle = " << cos_angle << std::endl;
			std::cout << "cos_angle on mesh = " << t10*t12*(t4*t5*t6*t20*(t25-l1*(t9+t23-t24-ep_f_v2_y*t8))+t4*t5*t6*t26*(l1*(t19+t29-t30-ep_f_v2_z*t8)-l2*t28))*(v1_x-v2_x)+t10*t12*(t4*t5*t6*t26*(t17-l1*(t3+t15-ep_f_v2_x*t8-ep_0_t*v1_x))-t4*t5*t6*t18*(t25-l1*(t9+t23-ep_f_v2_y*t8-ep_0_t*v1_y)))*(v1_z-v2_z)-t10*t12*(v1_y-v2_y)*(t4*t5*t6*t20*(t17-l1*(t3+t15-t16-ep_f_v2_x*t8))-t4*t5*t6*t18*(l2*t28-l1*(t19+t29-ep_f_v2_z*t8-ep_0_t*v1_z))) << std::endl;
			//std::cout << "normalized it's: " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z))/e_len << std::endl;
			//exit(1);
		};
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;

		int vnum = x.rows()/3;
		Eigen::VectorXd vals(fold_params.size()); 
		for (int i = 0; i < fold_params.size(); i++) {
			int v1_i(fold_params[i].v1),v2_i(fold_params[i].v2);
			int w1_i(fold_params[i].w1),w2_i(fold_params[i].w2);
			EdgePoint ep_b(fold_params[i].ep_b), ep_f(fold_params[i].ep_f);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));
			const double ep_0_t(fold_params[i].edge_t);
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double cos_angle = cos_angles[i];
			double fold_e_1 = folds_e1_l[i], fold_e_2 = folds_e2_l[i];
			double l1 = curve_l1[i]; double l2 = curve_l2[i];
			double fold_e_crease_angle = fold_e_crease_angles[i];

			double t2 = 1.0/fold_e_1;
			double t3 = 1.0/fold_e_2;
			double t4 = 1.0/l1;
			double t5 = sin(fold_e_crease_angle);
			double t6 = 1.0/t5;
			double t7 = w1_z-w2_z;
			double t8 = v1_z-v2_z;
			double t9 = v1_x-v2_x;
			double t10 = w1_y-w2_y;
			double t11 = v1_y-v2_y;
			double t12 = w1_x-w2_x;
			double t13 = ep_b_t-1.0;
			double t14 = 1.0/l2;
			double t15 = ep_f_t-1.0;
			double t16 = ep_0_t-1.0;
			double t17 = t16*v2_y;
			double t18 = t16*v2_z;
			double t19 = ep_0_t*l1;
			double t27 = ep_0_t*l2;
			double t20 = t19-t27;
			double t21 = t16*v2_x;
			double t22 = ep_b_t*ep_b_v1_z;
			double t25 = ep_0_t*v1_z;
			double t40 = ep_b_v2_z*t13;
			double t23 = t18+t22-t25-t40;
			double t24 = ep_f_t*ep_f_v1_z;
			double t26 = l2*t23;
			double t28 = ep_b_t*ep_b_v1_x;
			double t31 = ep_0_t*v1_x;
			double t45 = ep_b_v2_x*t13;
			double t29 = t21+t28-t31-t45;
			double t30 = ep_f_t*ep_f_v1_x;
			double t32 = l2*t29;
			double t33 = ep_b_t*ep_b_v1_y;
			double t36 = ep_0_t*v1_y;
			double t37 = ep_b_v2_y*t13;
			double t34 = t17+t33-t36-t37;
			double t35 = ep_f_t*ep_f_v1_y;
			double t38 = l2*t34;
			double t53 = ep_f_v2_y*t15;
			double t39 = t17+t35-t36-t53;
			double t48 = ep_f_v2_z*t15;
			double t41 = t18+t24-t25-t48;
			double t49 = l1*t41;
			double t42 = t26-t49;
			double t43 = l1*t16;
			double t50 = l2*t16;
			double t44 = t43-t50;
			double t51 = ep_f_v2_x*t15;
			double t46 = t21+t30-t31-t51;
			double t52 = l1*t46;
			double t47 = t32-t52;
			double t55 = l1*t39;
			double t54 = t38-t55;
			double t56 = t2*t3*t4*t6*t11*t14*t42;
			double t57 = t2*t3*t4*t6*t8*t14*t47;
			double t58 = t2*t3*t4*t6*t9*t14*t54;
			
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v1_i, -ep_b_t*t2*t3*t4*t6*t7*t11+ep_b_t*t2*t3*t4*t6*t8*t10);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v1_i+vnum, ep_b_t*t2*t3*t4*t6*t7*t9-ep_b_t*t2*t3*t4*t6*t8*t12);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v1_i+2*vnum, -ep_b_t*t2*t3*t4*t6*t9*t10+ep_b_t*t2*t3*t4*t6*t11*t12);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v2_i, t2*t3*t4*t6*t7*t11*t13-t2*t3*t4*t6*t8*t10*t13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v2_i+vnum, -t2*t3*t4*t6*t7*t9*t13+t2*t3*t4*t6*t8*t12*t13);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_b_v2_i+2*vnum, t2*t3*t4*t6*t9*t10*t13-t2*t3*t4*t6*t11*t12*t13);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v1_i, ep_f_t*t2*t3*t6*t7*t11*t14-ep_f_t*t2*t3*t6*t8*t10*t14);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v1_i+vnum, -ep_f_t*t2*t3*t6*t7*t9*t14+ep_f_t*t2*t3*t6*t8*t12*t14);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v1_i+2*vnum, ep_f_t*t2*t3*t6*t9*t10*t14-ep_f_t*t2*t3*t6*t11*t12*t14);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v2_i, -t2*t3*t6*t7*t11*t14*t15+t2*t3*t6*t8*t10*t14*t15);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v2_i+vnum, t2*t3*t6*t7*t9*t14*t15-t2*t3*t6*t8*t12*t14*t15);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, ep_f_v2_i+2*vnum, -t2*t3*t6*t9*t10*t14*t15+t2*t3*t6*t11*t12*t14*t15);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i, -t2*t6*(t3*t4*t10*t14*(t26-l1*(t18+t24-ep_f_v2_z*t15-ep_0_t*v1_z))-t3*t4*t7*t14*(t38-l1*(t17+t35-ep_f_v2_y*t15-ep_0_t*v1_y)))-t2*t3*t4*t6*t7*t11*t14*t20+t2*t3*t4*t6*t8*t10*t14*t20);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+vnum, -t2*t6*(t3*t4*t7*t14*(t32-l1*(t21+t30-ep_f_v2_x*t15-ep_0_t*v1_x))-t3*t4*t12*t14*t42)+t2*t3*t4*t6*t7*t9*t14*t20-t2*t3*t4*t6*t8*t12*t14*t20);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+2*vnum, -t2*t6*(t3*t4*t12*t14*(l2*t34-l1*t39)-t3*t4*t10*t14*t47)-t2*t3*t4*t6*t9*t10*t14*t20+t2*t3*t4*t6*t11*t12*t14*t20);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i, t2*t6*(t3*t4*t10*t14*t42-t3*t4*t7*t14*t54)+t2*t3*t4*t6*t7*t11*t14*t44-t2*t3*t4*t6*t8*t10*t14*t44);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+vnum, t2*t6*(t3*t4*t7*t14*t47-t3*t4*t12*t14*t42)-t2*t3*t4*t6*t7*t9*t14*t44+t2*t3*t4*t6*t8*t12*t14*t44);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+2*vnum, -t2*t6*(t3*t4*t10*t14*t47-t3*t4*t12*t14*t54)+t2*t3*t4*t6*t9*t10*t14*t44-t2*t3*t4*t6*t11*t12*t14*t44);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i, t56-t2*t3*t4*t6*t8*t14*t54);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i+vnum, t57-t2*t3*t4*t6*t9*t14*t42);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i+2*vnum, t58-t2*t3*t4*t6*t11*t14*t47);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i, -t56+t2*t3*t4*t6*t8*t14*t54);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i+vnum, -t57+t2*t3*t4*t6*t9*t14*t42);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i+2*vnum, -t58+t2*t3*t4*t6*t11*t14*t47);

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
		//std::cout << "const_cnt = " << const_cnt << std::endl;
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// The constraints are not linear but for now we use here gauss-newton..
	};
	const std::vector<MVTangentCreaseFold>& fold_params;
	std::vector<double> cos_angles;
	std::vector<double> folds_e1_l,folds_e2_l;
	std::vector<double> curve_l1,curve_l2;
	std::vector<double> fold_e_crease_angles;
};