#include "CurvedFoldingBiasObjective.h"

void CurvedFoldingBiasObjective::reset_folds() {
	curvedFoldBiases.clear();
	IJV.resize(0);
	is_H_cached = false;	
}

void CurvedFoldingBiasObjective::add_fold_bias(const CurvedFoldBias& foldBias) {
	curvedFoldBiases.push_back(foldBias);
	IJV.resize(curvedFoldBiases.size()*144);
	is_H_cached = false; // recompute H caching next time
}

double CurvedFoldingBiasObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	int vnum = x.rows()/3;

	int h_cnt = 0;
	for (auto curvedFold : curvedFoldBiases) {
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2), w1_i(curvedFold.w1), w2_i(curvedFold.w2);

		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
		const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
		const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

		const double ep_0_t(curvedFold.edge_t);

		double t2 = ep_0_t-1.0;
		double t3 = ep_b_t-1.0;
		double t4 = t2*v2_y;
		double t5 = ep_f_t-1.0;
		double t6 = t2*v2_x;
		double t7 = ep_b_t*ep_b_v1_x;
		double t11 = ep_0_t*v1_x;
		double t20 = ep_b_v2_x*t3;
		double t8 = t6+t7-t11-t20;
		double t9 = t2*v2_z;
		double t10 = ep_f_t*ep_f_v1_x;
		double t12 = ep_b_t*ep_b_v1_y;
		double t19 = ep_0_t*v1_y;
		double t22 = ep_b_v2_y*t3;
		double t13 = t4+t12-t19-t22;
		double t14 = ep_f_t*ep_f_v1_z;
		double t17 = ep_0_t*v1_z;
		double t24 = ep_f_v2_z*t5;
		double t15 = t9+t14-t17-t24;
		double t16 = ep_b_t*ep_b_v1_z;
		double t18 = ep_f_t*ep_f_v1_y;
		double t30 = ep_f_v2_y*t5;
		double t21 = t4+t18-t19-t30;
		double t27 = ep_f_v2_x*t5;
		double t23 = t6+t10-t11-t27;
		double t25 = t8*t15;
		double t29 = ep_b_v2_z*t3;
		double t26 = t9+t16-t17-t29;
		double t28 = t13*t15;
		double t31 = t28-t21*t26;
		double t32 = (t25-t23*t26)*(w1_y-w2_y)+(t25-t23*(t9+t16-ep_b_v2_z*t3-ep_0_t*v1_z))*(v1_y-v2_y)-t31*(v1_x-v2_x)-t31*(w1_x-w2_x)+(t13*(t6+t10-ep_f_v2_x*t5-ep_0_t*v1_x)-t8*(t4+t18-ep_f_v2_y*t5-ep_0_t*v1_y))*(v1_z-v2_z)-(w1_z-w2_z)*(t8*t21-t13*t23);
		e += t32*t32;


		if (dbg) {
			{
				/// TODO remove me
				double t2 = ep_0_t-1.0;
				double t3 = ep_b_t-1.0;
				double t4 = t2*v2_z;
				double t5 = ep_f_t-1.0;
				double t6 = t2*v2_y;
				double t7 = ep_f_t*ep_f_v1_z;
				double t10 = ep_0_t*v1_z;
				double t8 = t4+t7-t10-ep_f_v2_z*t5;
				double t9 = ep_b_t*ep_b_v1_z;
				double t11 = t2*v2_x;
				double t12 = ep_b_t*ep_b_v1_x;
				double t19 = ep_0_t*v1_x;
				double t13 = t11+t12-t19-ep_b_v2_x*t3;
				double t14 = ep_f_t*ep_f_v1_y;
				double t17 = ep_0_t*v1_y;
				double t15 = t6+t14-t17-ep_f_v2_y*t5;
				double t16 = ep_b_t*ep_b_v1_y;
				double t18 = ep_f_t*ep_f_v1_x;

				double B_fixed_x = -t8*(t6+t16-ep_b_v2_y*t3-ep_0_t*v1_y)+t15*(t4+t9-ep_b_v2_z*t3-ep_0_t*v1_z);
				double B_fixed_y = t8*t13-(t4+t9-t10-ep_b_v2_z*t3)*(t11+t18-ep_f_v2_x*t5-ep_0_t*v1_x);
				double B_fixed_z = (t6+t16-t17-ep_b_v2_y*t3)*(t11+t18-t19-ep_f_v2_x*t5)-t13*t15;

				Eigen::RowVector3d B; B<< B_fixed_x,B_fixed_y,B_fixed_z;
				std::cout << "B = " << B << std::endl;
				Eigen::RowVector3d e1; e1 << v1_x-v2_x,v1_y-v2_y,v1_z-v2_z;
				Eigen::RowVector3d e2; e2 << w1_x-w2_x,w1_y-w2_y,w1_z-w2_z;
				std::cout << "e1.dot(B) = " << e1.dot(B) << std::endl;
				std::cout << "e2.dot(B) = " << e2.dot(B) << std::endl;

				std::cout << "pow(e1.dot(B)+e2.dot(B),2) = " << pow(e1.dot(B)+e2.dot(B),2) << " but t32*t32 = " << t32*t32 << std::endl;
			}
		}
		

  }
  return e;
}

Eigen::VectorXd CurvedFoldingBiasObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;

  for (auto curvedFold : curvedFoldBiases) {
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2), w1_i(curvedFold.w1), w2_i(curvedFold.w2);

		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
		const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
		const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

		const double ep_0_t(curvedFold.edge_t);

		double t2 = ep_f_t-1.0;
		double t3 = ep_0_t-1.0;
		double t4 = ep_f_t*ep_f_v1_y;
		double t5 = t3*v2_y;
		double t11 = ep_0_t*v1_y;
		double t12 = ep_f_v2_y*t2;
		double t6 = t4+t5-t11-t12;
		double t7 = ep_f_t*ep_f_v1_z;
		double t8 = t3*v2_z;
		double t18 = ep_0_t*v1_z;
		double t19 = ep_f_v2_z*t2;
		double t9 = t7+t8-t18-t19;
		double t10 = v1_z-v2_z;
		double t13 = ep_b_t-1.0;
		double t14 = t3*v2_x;
		double t15 = v1_y-v2_y;
		double t16 = ep_b_t*ep_b_v1_x;
		double t21 = ep_0_t*v1_x;
		double t27 = ep_b_v2_x*t13;
		double t17 = t14+t16-t21-t27;
		double t20 = ep_f_t*ep_f_v1_x;
		double t22 = ep_b_t*ep_b_v1_y;
		double t29 = ep_b_v2_y*t13;
		double t23 = t5-t11+t22-t29;
		double t24 = ep_b_t*ep_b_v1_z;
		double t33 = ep_b_v2_z*t13;
		double t25 = t8-t18+t24-t33;
		double t26 = w1_z-w2_z;
		double t28 = t6*t17;
		double t34 = ep_f_v2_x*t2;
		double t30 = t14+t20-t21-t34;
		double t31 = w1_y-w2_y;
		double t32 = t9*t17;
		double t42 = t25*t30;
		double t35 = t32-t42;
		double t36 = t6*t25;
		double t44 = t9*t23;
		double t37 = t36-t44;
		double t41 = t23*t30;
		double t38 = t28-t41;
		double t39 = v1_x-v2_x;
		double t40 = w1_x-w2_x;
		double t43 = t15*t35;
		double t45 = t37*t39;
		double t46 = t31*t35;
		double t47 = t37*t40;
		double t49 = t10*t38;
		double t50 = t26*t38;
		double t48 = t43+t45+t46+t47-t49-t50;
		double t51 = ep_0_t*t6;
		double t57 = ep_0_t*t23;
		double t52 = t51-t57;
		double t53 = ep_0_t*t9;
		double t54 = ep_0_t*t17;
		double t56 = ep_0_t*t30;
		double t55 = t54-t56;
		double t58 = t3*t6;
		double t66 = t3*t23;
		double t59 = t58-t66;
		double t60 = t3*t9;
		double t62 = t3*t25;
		double t61 = t60-t62;
		double t63 = t3*t17;
		double t65 = t3*t30;
		double t64 = t63-t65;
		double t67 = -t36+t44;
		double t68 = t39*t67;
		double t69 = t40*t67;
		double t70 = -t43-t46+t49+t50+t68+t69;
		double t71 = t67*t70*2.0;
		double t72 = t38*t70*2.0;

		grad(ep_b_v1_i) += t48*(ep_b_t*t6*t10-ep_b_t*t9*t15+ep_b_t*t6*t26-ep_b_t*t9*t31)*-2.0;
		grad(ep_b_v1_i+vnum) += t48*(ep_b_t*t10*t30-ep_b_t*t9*t39-ep_b_t*t9*t40+ep_b_t*t26*t30)*2.0;
		grad(ep_b_v1_i+2*vnum) += t48*(ep_b_t*t6*t39-ep_b_t*t15*t30+ep_b_t*t6*t40-ep_b_t*t30*t31)*2.0;
		grad(ep_b_v2_i) += t48*(t6*t10*t13-t9*t13*t15+t6*t13*t26-t9*t13*t31)*2.0;
		grad(ep_b_v2_i+vnum) += t48*(t10*t13*t30-t9*t13*t39-t9*t13*t40+t13*t26*t30)*-2.0;
		grad(ep_b_v2_i+2*vnum) += t48*(t6*t13*t39-t13*t15*t30+t6*t13*t40-t13*t30*t31)*-2.0;
		grad(ep_f_v1_i) += (ep_f_t*t10*t23-ep_f_t*t15*t25+ep_f_t*t23*t26-ep_f_t*t25*t31)*(t43+t45+t46+t47-t49-t50)*2.0;
		grad(ep_f_v1_i+vnum) += t48*(ep_f_t*t10*t17+ep_f_t*t17*t26-ep_f_t*t25*t39-ep_f_t*t25*t40)*-2.0;
		grad(ep_f_v1_i+2*vnum) += (ep_f_t*t15*t17+ep_f_t*t17*t31-ep_f_t*t23*t39-ep_f_t*t23*t40)*(t43+t45+t46+t47-t49-t50)*2.0;
		grad(ep_f_v2_i) += t48*(t2*t10*t23-t2*t15*t25+t2*t23*t26-t2*t25*t31)*-2.0;
		grad(ep_f_v2_i+vnum) += (t2*t10*t17+t2*t17*t26-t2*t25*t39-t2*t25*t40)*(t43+t45+t46+t47-t49-t50)*2.0;
		grad(ep_f_v2_i+2*vnum) += t48*(t2*t15*t17+t2*t17*t31-t2*t23*t39-t2*t23*t40)*-2.0;
		grad(v1_i) += (t43+t45+t46+t47-t49-t50)*(t36-t44+t10*t52+t26*t52-t15*(t53-ep_0_t*t25)-t31*(t53-ep_0_t*t25))*2.0;
		grad(v1_i+vnum) += (t43+t45+t46+t47-t49-t50)*(t32-t42+t10*t55+t26*t55+t39*(t53-ep_0_t*t25)+t40*(t53-ep_0_t*t25))*2.0;
		grad(v1_i+2*vnum) += t48*(t28-t41+t15*t55+t31*t55+t39*t52+t40*t52)*-2.0;
		grad(v2_i) += t48*(t36-t44+t10*t59-t15*t61+t26*t59-t31*t61)*-2.0;
		grad(v2_i+vnum) += t48*(t32-t42+t10*t64+t26*t64+t39*t61+t40*t61)*-2.0;
		grad(v2_i+2*vnum) += (t28-t41+t15*t64+t31*t64+t39*(t58-t66)+t40*(t58-t66))*(t43+t45+t46+t47-t49-t50)*2.0;
		grad(w1_i) += t71;
		grad(w1_i+vnum) += t35*t70*-2.0;
		grad(w1_i+2*vnum) += t72;
		grad(w2_i) += -t71;
		grad(w2_i+vnum) += t35*t70*2.0;
		grad(w2_i+2*vnum) += -t72;

  }
  return grad;
}

void CurvedFoldingBiasObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0;

  int ijv_cnt = 0;
  #pragma clang loop vectorize(enable)
	 for (auto curvedFold : curvedFoldBiases) {
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2), w1_i(curvedFold.w1), w2_i(curvedFold.w2);

		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
		const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
		const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

		const double ep_0_t(curvedFold.edge_t);

		// Get the fixed binormal
		double B_fixed_x,B_fixed_y,B_fixed_z;
		{
			double t2 = ep_0_t-1.0;
			double t3 = ep_b_t-1.0;
			double t4 = t2*v2_z;
			double t5 = ep_f_t-1.0;
			double t6 = t2*v2_y;
			double t7 = ep_f_t*ep_f_v1_z;
			double t10 = ep_0_t*v1_z;
			double t8 = t4+t7-t10-ep_f_v2_z*t5;
			double t9 = ep_b_t*ep_b_v1_z;
			double t11 = t2*v2_x;
			double t12 = ep_b_t*ep_b_v1_x;
			double t19 = ep_0_t*v1_x;
			double t13 = t11+t12-t19-ep_b_v2_x*t3;
			double t14 = ep_f_t*ep_f_v1_y;
			double t17 = ep_0_t*v1_y;
			double t15 = t6+t14-t17-ep_f_v2_y*t5;
			double t16 = ep_b_t*ep_b_v1_y;
			double t18 = ep_f_t*ep_f_v1_x;

			B_fixed_x = -t8*(t6+t16-ep_b_v2_y*t3-ep_0_t*v1_y)+t15*(t4+t9-ep_b_v2_z*t3-ep_0_t*v1_z);
			B_fixed_y = t8*t13-(t4+t9-t10-ep_b_v2_z*t3)*(t11+t18-ep_f_v2_x*t5-ep_0_t*v1_x);
			B_fixed_z = (t6+t16-t17-ep_b_v2_y*t3)*(t11+t18-t19-ep_f_v2_x*t5)-t13*t15;
		}

		double t2 = B_fixed_x*B_fixed_x;
		double t3 = t2*2.0;
		double t4 = B_fixed_x*B_fixed_y*2.0;
		double t5 = B_fixed_x*B_fixed_z*2.0;
		double t6 = B_fixed_y*B_fixed_y;
		double t7 = t6*2.0;
		double t8 = B_fixed_y*B_fixed_z*2.0;
		double t9 = B_fixed_z*B_fixed_z;
		double t10 = t9*2.0;

		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i, t3));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i, t3));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i, t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i, t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i, t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i, t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i, t3));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i, t3));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i, t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i, t4));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i, t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i, t5));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i, t3));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i, t3));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i, t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i, t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i, t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i, t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i, t3));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i, -t3));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i+vnum, -t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i+2*vnum, -t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i, t3));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i+vnum, t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i+2*vnum, t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i, t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i, -t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i+vnum, -t7));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i+2*vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i, t4));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i+vnum, t7));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i+2*vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i, t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i+2*vnum, t10));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i, -t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i+vnum, -t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i+2*vnum, -t10));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i, t5));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i+vnum, t8));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i+2*vnum, t10));
  }
}