#include "CurvedFoldingBiasObjective.h"

void CurvedFoldingBiasObjective::reset_folds() {
	curvedFoldBiases.clear();
	IJV.resize(0);
	is_H_cached = false;	
}

void CurvedFoldingBiasObjective::add_fold_bias(const CurvedFoldBias& foldBias) {
	curvedFoldBiases.push_back(foldBias);
	IJV.resize(curvedFoldBiases.size()*36);
	is_H_cached = false; // recompute H caching next time
}

double CurvedFoldingBiasObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	int vnum = x.rows()/3;

	int h_cnt = 0;
	for (auto curvedFold : curvedFoldBiases) {
		int ep_0_v1_i(curvedFold.ep_0.edge.v1), ep_0_v2_i(curvedFold.ep_0.edge.v2); const double ep_0_t(curvedFold.ep_0.t);
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2);

		const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
		const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));

		double t2 = ep_0_t-1.0;
		double t3 = ep_0_t*ep_0_v1_y;
		double t4 = ep_b_t-1.0;
		double t5 = ep_0_t*ep_0_v1_x;
		double t6 = ep_f_t-1.0;
		double t7 = ep_b_v2_x*t4;
		double t13 = ep_0_v2_x*t2;
		double t17 = ep_b_t*ep_b_v1_x;
		double t8 = t5+t7-t13-t17;
		double t9 = ep_f_v2_y*t6;
		double t11 = ep_0_v2_y*t2;
		double t26 = ep_f_t*ep_f_v1_y;
		double t10 = t3+t9-t11-t26;
		double t12 = ep_b_v2_y*t4;
		double t14 = ep_f_v2_x*t6;
		double t15 = ep_0_v2_z*t2;
		double t16 = t8*t10;
		double t23 = ep_f_t*ep_f_v1_x;
		double t18 = t5-t13+t14-t23;
		double t19 = ep_f_t*ep_f_v1_z;
		double t21 = ep_0_t*ep_0_v1_z;
		double t25 = ep_f_v2_z*t6;
		double t20 = t15+t19-t21-t25;
		double t22 = ep_b_t*ep_b_v1_z;
		double t27 = ep_b_t*ep_b_v1_y;
		double t24 = t3-t11+t12-t27;
		double t28 = t10*(t15-t21+t22-ep_b_v2_z*t4);
		double t29 = t28-t20*t24;
		double t30 = (t18*(t15-t21+t22-ep_b_v2_z*t4)-t8*t20)*(-t3+t11+v2_y)+(t16-(t3+t12-ep_b_t*ep_b_v1_y-ep_0_v2_y*t2)*(t5+t14-ep_f_t*ep_f_v1_x-ep_0_v2_x*t2))*(t15+v1_z-ep_0_t*ep_0_v1_z)+(t8*t20-t18*(t15+t22-ep_0_t*ep_0_v1_z-ep_b_v2_z*t4))*(-t3+t11+v1_y)-(t16-t18*t24)*(t15+v2_z-ep_0_t*ep_0_v1_z)+t29*(-t5+t13+v1_x)-t29*(-t5+t13+v2_x);
		e += t30*t30;

  }
  return e;
}

Eigen::VectorXd CurvedFoldingBiasObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;

  for (auto curvedFold : curvedFoldBiases) {
		int ep_0_v1_i(curvedFold.ep_0.edge.v1), ep_0_v2_i(curvedFold.ep_0.edge.v2); const double ep_0_t(curvedFold.ep_0.t);
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2);

		const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
		const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));

		double t2 = ep_0_t*ep_0_v1_y;
		double t3 = ep_0_t-1.0;
		double t4 = ep_b_t-1.0;
		double t5 = ep_b_v2_y*t4;
		double t8 = ep_0_v2_y*t3;
		double t17 = ep_b_t*ep_b_v1_y;
		double t6 = t2+t5-t8-t17;
		double t7 = ep_0_t*t6;
		double t9 = ep_f_t-1.0;
		double t10 = ep_f_v2_y*t9;
		double t11 = ep_0_v2_z*t3;
		double t12 = ep_b_t*ep_b_v1_z;
		double t14 = ep_0_t*ep_0_v1_z;
		double t26 = ep_b_v2_z*t4;
		double t13 = t11+t12-t14-t26;
		double t15 = ep_f_t*ep_f_v1_z;
		double t21 = ep_f_t*ep_f_v1_y;
		double t16 = t2-t8+t10-t21;
		double t18 = ep_0_t*ep_0_v1_x;
		double t19 = ep_b_v2_x*t4;
		double t22 = ep_0_v2_x*t3;
		double t24 = ep_b_t*ep_b_v1_x;
		double t20 = t18+t19-t22-t24;
		double t23 = ep_f_v2_x*t9;
		double t29 = ep_f_v2_z*t9;
		double t25 = t11-t14+t15-t29;
		double t31 = ep_f_t*ep_f_v1_x;
		double t27 = t18-t22+t23-t31;
		double t28 = -t2+t8+v1_y;
		double t30 = t20*t25;
		double t32 = -t2+t8+v2_y;
		double t33 = ep_0_t*t13;
		double t47 = ep_0_t*t25;
		double t34 = t33-t47;
		double t56 = ep_0_t*t16;
		double t35 = t7-t56;
		double t36 = t11-t14+v2_z;
		double t37 = t6*t27;
		double t38 = t11-t14+v1_z;
		double t48 = t16*t20;
		double t39 = t37-t48;
		double t40 = t13*t16;
		double t50 = t6*t25;
		double t41 = t40-t50;
		double t49 = t13*t27;
		double t42 = t30-t49;
		double t43 = ep_0_t*t20;
		double t55 = ep_0_t*t27;
		double t44 = t43-t55;
		double t45 = -t18+t22+v1_x;
		double t46 = -t18+t22+v2_x;
		double t51 = t41*t45;
		double t52 = t28*t42;
		double t53 = t36*t39;
		double t57 = t38*t39;
		double t58 = t32*t42;
		double t59 = t41*t46;
		double t54 = t51+t52+t53-t57-t58-t59;
		double t60 = t3*t6;
		double t68 = t3*t16;
		double t61 = t60-t68;
		double t62 = t3*t13;
		double t66 = t3*t25;
		double t63 = t62-t66;
		double t64 = t3*t20;
		double t67 = t3*t27;
		double t65 = t64-t67;
		double t69 = t36*(t37-t48);
		double t71 = t27*(t11+t12-t14-t26);
		double t72 = t30-t71;
		double t70 = t28*t72;
		double t74 = t16*(t11+t12-t14-t26);
		double t75 = t50-t74;
		double t73 = t46*t75;
		double t76 = t32*t72;
		double t77 = t45*t75;
		double t78 = t57-t69-t70-t73+t76+t77;
		double t79 = -t50+t74;
		double t80 = -t30+t71;
		double t81 = t46*t79;
		double t82 = -t37+t48;
		double t83 = t28*t80;
		double t84 = t36*t82;
		double t86 = t38*t82;
		double t87 = t32*t80;
		double t88 = t45*t79;
		double t85 = t81+t83+t84-t86-t87-t88;

		grad(ep_0_v1_i) += t54*(t28*t34-t32*t34-t35*t36+t35*t38)*-2.0;
		grad(ep_0_v1_i+v_num) += t54*(t34*t45-t34*t46-t36*t44+t38*t44)*2.0;
		grad(ep_0_v1_i+2*v_num) += t54*(t28*t44-t32*t44-t35*t45+t35*t46)*-2.0;
		grad(ep_0_v2_i) += t54*(t28*t63-t32*t63-t36*t61+t38*t61)*2.0;
		grad(ep_0_v2_i+v_num) += (t36*t65-t38*t65+t46*t63-t45*(t62-t3*(t11-t14+t15-t29)))*(t51-t57-t58+t69+t70+t73)*2.0;
		grad(ep_0_v2_i+2*vnum) += (t28*t65-t32*t65-t45*t61+t46*t61)*(t51-t57+t69+t70+t73-t32*t72)*2.0;
		grad(ep_b_v1_i) += (ep_b_t*t28*(t11-t14+t15-t29)-ep_b_t*t16*t36+ep_b_t*t16*t38-ep_b_t*t25*t32)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(ep_b_v1_i+v_num) += t78*(ep_b_t*t27*t36-ep_b_t*t27*t38-ep_b_t*t25*t45+ep_b_t*t25*t46)*2.0;
		grad(ep_b_v1_i+2*v_num) += (ep_b_t*t27*t28-ep_b_t*t27*t32-ep_b_t*t16*t45+ep_b_t*t16*t46)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(ep_b_v2_i) += t78*(t4*t16*t36-t4*t25*t28-t4*t16*t38+t4*t25*t32)*2.0;
		grad(ep_b_v2_i+v_num) += (t4*t45*(t11-t14+t15-t29)-t4*t27*t36+t4*t27*t38-t4*t25*t46)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(ep_b_v2_i+2*vnum) += t78*(t4*t27*t28-t4*t27*t32-t4*t16*t45+t4*t16*t46)*-2.0;
		grad(ep_f_v1_i) += t78*(ep_f_t*t13*t28-ep_f_t*t6*t36+ep_f_t*t6*t38-ep_f_t*t13*t32)*-2.0;
		grad(ep_f_v1_i+v_num) += (ep_f_t*t45*(t11+t12-t14-t26)-ep_f_t*t20*t36+ep_f_t*t20*t38-ep_f_t*t13*t46)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(ep_f_v1_i+2*v_num) += t78*(ep_f_t*t20*t28-ep_f_t*t6*t45+ep_f_t*t6*t46-ep_f_t*t20*t32)*-2.0;
		grad(ep_f_v2_i) += (t9*t28*(t11+t12-t14-t26)-t6*t9*t36+t6*t9*t38-t9*t13*t32)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(ep_f_v2_i+v_num) += t78*(t9*t20*t36-t9*t13*t45-t9*t20*t38+t9*t13*t46)*2.0;
		grad(ep_f_v2_i+2*v_num) += (t9*t20*t28-t6*t9*t45+t6*t9*t46-t9*t20*t32)*(t57-t69+t76+t77-t28*(t30-t71)-t46*(t50-t74))*2.0;
		grad(v1_i) += t79*(t57-t69-t70+t76+t81-t45*t79)*-2.0;
		grad(v1_i+v_num) += t80*(t57-t69+t81+t83-t32*t80-t45*t79)*2.0;
		grad(v1_i+2*v_num) += t82*t85*-2.0;
		grad(v2_i) += t79*t85*2.0;
		grad(v2_i+v_num) += t80*t85*-2.0;
		grad(v2_i+2*v_num) += t82*(t81+t83+t84-t86-t87-t88)*2.0;
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
		int ep_0_v1_i(curvedFold.ep_0.edge.v1), ep_0_v2_i(curvedFold.ep_0.edge.v2); const double ep_0_t(curvedFold.ep_0.t);
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2);

		const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
		const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
		const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
		const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
		const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
		const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

		const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
		const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));

		// Get the fixed binormal
		double B_fixed_x,B_fixed_y,B_fixed_z;
		{
			double t2 = ep_0_t-1.0;
			double t3 = ep_0_t*ep_0_v1_z;
			double t4 = ep_b_t-1.0;
			double t5 = ep_0_t*ep_0_v1_y;
			double t6 = ep_f_t-1.0;
			double t7 = ep_f_v2_z*t6;
			double t9 = ep_0_v2_z*t2;
			double t8 = t3+t7-t9-ep_f_t*ep_f_v1_z;
			double t10 = ep_b_v2_z*t4;
			double t11 = ep_0_t*ep_0_v1_x;
			double t12 = ep_b_v2_x*t4;
			double t18 = ep_0_v2_x*t2;
			double t13 = t11+t12-t18-ep_b_t*ep_b_v1_x;
			double t14 = ep_f_v2_y*t6;
			double t16 = ep_0_v2_y*t2;
			double t15 = t5+t14-t16-ep_f_t*ep_f_v1_y;
			double t17 = ep_b_v2_y*t4;
			double t19 = ep_f_v2_x*t6;

			B_fixed_x = -t8*(t5+t17-ep_b_t*ep_b_v1_y-ep_0_v2_y*t2)+t15*(t3+t10-ep_b_t*ep_b_v1_z-ep_0_v2_z*t2);
  			B_fixed_y = t8*t13-(t3-t9+t10-ep_b_t*ep_b_v1_z)*(t11+t19-ep_f_t*ep_f_v1_x-ep_0_v2_x*t2);
  			B_fixed_z = -t13*t15+(t5-t16+t17-ep_b_t*ep_b_v1_y)*(t11-t18+t19-ep_f_t*ep_f_v1_x);
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

		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i, t3);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i+v_num, t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i+2*v_num, t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i, -t3);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i+v_num, -t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i+2*v_num, -t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i, t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i+v_num, t7);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i+2*v_num, t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i, -t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i+v_num, -t7);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i+2*v_num, -t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i, t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i+v_num, t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i+2*v_num, t10);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i, -t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i+v_num, -t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i+2*v_num, -t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i, -t3);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i+v_num, -t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i+2*v_num, -t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i, t3);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i+v_num, t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i+2*v_num, t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i, -t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i+v_num, -t7);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i+2*v_num, -t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i, t4);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i+v_num, t7);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i+2*v_num, t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i, -t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+v_num, -t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+2*v_num, -t10);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v2_i, t5);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v2_i+v_num, t8);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+2*v_num, t10);

      h_cnt++;
  }
}

/*
Hessian
A0[ v1_x][ v1_x] = t3; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i, t3);
  A0[ v1_x][ v1_y] = t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i+v_num, t4);
  A0[ v1_x][ v1_z] = t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v1_i+2*v_num, t5);
  A0[ v1_x][ v2_x] = -t3; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i, -t3);
  A0[ v1_x][ v2_y] = -t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i+v_num, -t4);
  A0[ v1_x][ v2_z] = -t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i,v2_i+2*v_num, -t5);
  A0[ v1_y][ v1_x] = t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i, t4);
  A0[ v1_y][ v1_y] = t7; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i+v_num, t7);
  A0[ v1_y][ v1_z] = t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v1_i+2*v_num, t8);
  A0[ v1_y][ v2_x] = -t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i, -t4);
  A0[ v1_y][ v2_y] = -t7; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i+v_num, -t7);
  A0[ v1_y][ v2_z] = -t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+v_num,v2_i+2*v_num, -t8);
  A0[ v1_z][ v1_x] = t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i, t5);
  A0[ v1_z][ v1_y] = t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i+v_num, t8);
  A0[ v1_z][ v1_z] = t10; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v1_i+2*v_num, t10);
  A0[ v1_z][ v2_x] = -t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i, -t5);
  A0[ v1_z][ v2_y] = -t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i+v_num, -t8);
  A0[ v1_z][ v2_z] = -t10; IJV[ijv_cnt++] = Eigen::Triplet<double>(v1_i+2*v_num,v2_i+2*v_num, -t5);
  A0[ v2_x][ v1_x] = -t3; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i, -t3);
  A0[ v2_x][ v1_y] = -t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i+v_num, -t4);
  A0[ v2_x][ v1_z] = -t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v1_i+2*v_num, -t5);
  A0[ v2_x][ v2_x] = t3; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i, t3);
  A0[ v2_x][ v2_y] = t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i+v_num, t4);
  A0[ v2_x][ v2_z] = t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i,v2_i+2*v_num, t5);
  A0[ v2_y][ v1_x] = -t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i, -t4);
  A0[ v2_y][ v1_y] = -t7; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i+v_num, -t7);
  A0[ v2_y][ v1_z] = -t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v1_i+2*v_num, -t8);
  A0[ v2_y][ v2_x] = t4; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i, t4);
  A0[ v2_y][ v2_y] = t7; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i+v_num, t7);
  A0[ v2_y][ v2_z] = t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+v_num,v2_i+2*v_num, t8);
  A0[ v2_z][ v1_x] = -t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i, -t5);
  A0[ v2_z][ v1_y] = -t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+v_num, -t8);
  A0[ v2_z][ v1_z] = -t10; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+2*v_num, -t10);
  A0[ v2_z][ v2_x] = t5; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v2_i, t5);
  A0[ v2_z][ v2_y] = t8; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v2_i+v_num, t8);
  A0[ v2_z][ v2_z] = t10; IJV[ijv_cnt++] = Eigen::Triplet<double>(v2_i+2*v_num,v1_i+2*v_num, t10);
*/