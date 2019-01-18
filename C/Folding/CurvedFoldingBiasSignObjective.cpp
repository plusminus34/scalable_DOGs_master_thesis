#include "CurvedFoldingBiasSignObjective.h"

void CurvedFoldingBiasSignObjective::reset_folds() {
	curvedFoldBiases.clear();
	IJV.resize(0);
	is_H_cached = false;	
}

void CurvedFoldingBiasSignObjective::add_fold_bias(const CurvedFoldBias& foldBias) {
	curvedFoldBiases.push_back(foldBias);
	IJV.resize(curvedFoldBiases.size()*144);
	is_H_cached = false; // recompute H caching next time
}

double CurvedFoldingBiasSignObjective::obj(const Eigen::VectorXd& x) const {
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
		double t23 = ep_b_v2_x*t3;
		double t8 = t6+t7-t11-t23;
		double t9 = t2*v2_z;
		double t10 = ep_f_t*ep_f_v1_x;
		double t12 = ep_b_t*ep_b_v1_y;
		double t19 = ep_0_t*v1_y;
		double t26 = ep_b_v2_y*t3;
		double t13 = t4+t12-t19-t26;
		double t14 = ep_f_t*ep_f_v1_z;
		double t17 = ep_0_t*v1_z;
		double t29 = ep_f_v2_z*t5;
		double t15 = t9+t14-t17-t29;
		double t16 = ep_b_t*ep_b_v1_z;
		double t18 = ep_f_t*ep_f_v1_y;
		double t27 = ep_f_v2_x*t5;
		double t20 = t6+t10-t11-t27;
		double t31 = ep_b_v2_z*t3;
		double t21 = t9+t16-t17-t31;
		double t24 = ep_f_v2_y*t5;
		double t22 = t4+t18-t19-t24;
		double t25 = t8*t22;
		double t28 = t25-t13*t20;
		double t30 = t8*t15;
		double t32 = t30-t20*t21;
		double t33 = t13*t15;
		double t34 = t33-t21*t22;
		double t35 = tanh(alpha*(t28*(v1_z-v2_z)-t32*(v1_y-v2_y)+t34*(v1_x-v2_x)))+tanh(alpha*(t28*(w1_z-w2_z)-t32*(w1_y-w2_y)+t34*(w1_x-w2_x)));
		e += t35*t35;


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

				std::cout << "pow(e1.dot(B)+e2.dot(B),2) = " << pow(e1.dot(B)+e2.dot(B),2) << " but e = " << t35*t35 << std::endl;
			}
		}
		

  }
  return e;
}

Eigen::VectorXd CurvedFoldingBiasSignObjective::grad(const Eigen::VectorXd& x) const {
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

		double t2 = ep_0_t-1.0;
		double t3 = ep_b_t-1.0;
		double t4 = t2*v2_y;
		double t5 = ep_f_t-1.0;
		double t6 = t2*v2_x;
		double t7 = ep_b_t*ep_b_v1_x;
		double t11 = ep_0_t*v1_x;
		double t23 = ep_b_v2_x*t3;
		double t8 = t6+t7-t11-t23;
		double t9 = t2*v2_z;
		double t10 = ep_f_t*ep_f_v1_x;
		double t12 = ep_b_t*ep_b_v1_y;
		double t19 = ep_0_t*v1_y;
		double t26 = ep_b_v2_y*t3;
		double t13 = t4+t12-t19-t26;
		double t14 = ep_f_t*ep_f_v1_z;
		double t17 = ep_0_t*v1_z;
		double t29 = ep_f_v2_z*t5;
		double t15 = t9+t14-t17-t29;
		double t16 = ep_b_t*ep_b_v1_z;
		double t18 = ep_f_t*ep_f_v1_y;
		double t27 = ep_f_v2_x*t5;
		double t20 = t6+t10-t11-t27;
		double t31 = ep_b_v2_z*t3;
		double t21 = t9+t16-t17-t31;
		double t24 = ep_f_v2_y*t5;
		double t22 = t4+t18-t19-t24;
		double t25 = t8*t22;
		double t37 = t13*t20;
		double t28 = t25-t37;
		double t30 = t8*t15;
		double t39 = t20*t21;
		double t32 = t30-t39;
		double t33 = t13*t15;
		double t41 = t21*t22;
		double t34 = t33-t41;
		double t35 = v1_z-v2_z;
		double t36 = v1_y-v2_y;
		double t38 = t28*t35;
		double t40 = v1_x-v2_x;
		double t42 = t34*t40;
		double t54 = t32*t36;
		double t43 = t38+t42-t54;
		double t44 = alpha*t43;
		double t45 = tanh(t44);
		double t46 = w1_z-w2_z;
		double t47 = w1_y-w2_y;
		double t48 = t28*t46;
		double t49 = w1_x-w2_x;
		double t50 = t34*t49;
		double t55 = t32*t47;
		double t51 = t48+t50-t55;
		double t52 = alpha*t51;
		double t53 = tanh(t52);
		double t56 = t45+t53;
		double t57 = t45*t45;
		double t58 = t57-1.0;
		double t59 = t53*t53;
		double t60 = t59-1.0;
		double t61 = ep_0_t*t13;
		double t69 = ep_0_t*t22;
		double t62 = t61-t69;
		double t63 = ep_0_t*t15;
		double t65 = ep_0_t*t21;
		double t64 = t63-t65;
		double t66 = ep_0_t*t8;
		double t68 = ep_0_t*t20;
		double t67 = t66-t68;
		double t70 = t2*t13;
		double t78 = t2*t22;
		double t71 = t70-t78;
		double t72 = t2*t15;
		double t73 = t2*t8;
		double t77 = t2*t20;
		double t74 = t73-t77;
		double t76 = t2*t21;
		double t75 = t72-t76;
		double t79 = alpha*t32*t56*t60*2.0;

		grad(ep_b_v1_i) += t56*(alpha*t58*(ep_b_t*t15*t36-ep_b_t*t22*t35)+alpha*t60*(ep_b_t*t15*t47-ep_b_t*t22*t46))*2.0;
		grad(ep_b_v1_i+vnum) += t56*(alpha*t58*(ep_b_t*t15*t40-ep_b_t*t20*t35)+alpha*t60*(ep_b_t*t15*t49-ep_b_t*t20*t46))*-2.0;
		grad(ep_b_v1_i+2*vnum) += t56*(alpha*t58*(ep_b_t*t20*t36-ep_b_t*t22*t40)+alpha*t60*(ep_b_t*t20*t47-ep_b_t*t22*t49))*-2.0;
		grad(ep_b_v2_i) += t56*(alpha*t58*(t3*t15*t36-t3*t22*t35)+alpha*t60*(t3*t15*t47-t3*t22*t46))*-2.0;
		grad(ep_b_v2_i+vnum) += t56*(alpha*t58*(t3*t15*t40-t3*t20*t35)+alpha*t60*(t3*t15*t49-t3*t20*t46))*2.0;
		grad(ep_b_v2_i+2*vnum) += t56*(alpha*t58*(t3*t20*t36-t3*t22*t40)+alpha*t60*(t3*t20*t47-t3*t22*t49))*2.0;
		grad(ep_f_v1_i) += t56*(alpha*t58*(ep_f_t*t13*t35-ep_f_t*t21*t36)+alpha*t60*(ep_f_t*t13*t46-ep_f_t*t21*t47))*2.0;
		grad(ep_f_v1_i+vnum) += t56*(alpha*t58*(ep_f_t*t8*t35-ep_f_t*t21*t40)+alpha*t60*(ep_f_t*t8*t46-ep_f_t*t21*t49))*-2.0;
		grad(ep_f_v1_i+2*vnum) += t56*(alpha*t58*(ep_f_t*t8*t36-ep_f_t*t13*t40)+alpha*t60*(ep_f_t*t8*t47-ep_f_t*t13*t49))*2.0;
		grad(ep_f_v2_i) += t56*(alpha*t58*(t5*t13*t35-t5*t21*t36)+alpha*t60*(t5*t13*t46-t5*t21*t47))*-2.0;
		grad(ep_f_v2_i+vnum) += t56*(alpha*t58*(t5*t8*t35-t5*t21*t40)+alpha*t60*(t5*t8*t46-t5*t21*t49))*2.0;
		grad(ep_f_v2_i+2*vnum) += t56*(alpha*t58*(t5*t8*t36-t5*t13*t40)+alpha*t60*(t5*t8*t47-t5*t13*t49))*-2.0;
		grad(v1_i) += t56*(alpha*t60*(t46*t62+t47*t64)+alpha*t58*(t33-t41+t35*t62+t36*t64))*-2.0;
		grad(v1_i+vnum) += t56*(alpha*t60*(t46*t67+t49*(t63-t65))+alpha*t58*(t30-t39+t35*t67+t40*(t63-t65)))*2.0;
		grad(v1_i+2*vnum) += t56*(alpha*t60*(t49*t62-t47*t67)-alpha*t58*(t25-t37-t40*t62+t36*t67))*2.0;
		grad(v2_i) += t56*(alpha*t60*(t46*t71+t47*t75)+alpha*t58*(t33-t41+t35*t71+t36*t75))*2.0;
		grad(v2_i+vnum) += t56*(alpha*t60*(t46*t74+t49*t75)+alpha*t58*(t30-t39+t35*t74+t40*t75))*-2.0;
		grad(v2_i+2*vnum) += t56*(alpha*t60*(t49*t71-t47*t74)-alpha*t58*(t25-t37+t36*t74-t40*t71))*-2.0;
		grad(w1_i) += alpha*t34*t56*t60*-2.0;
		grad(w1_i+vnum) += t79;
		grad(w1_i+2*vnum) += alpha*t28*t56*t60*-2.0;
		grad(w2_i) += alpha*t34*t56*t60*2.0;
		grad(w2_i+vnum) += -t79;
		grad(w2_i+2*vnum) += alpha*t28*t56*t60*2.0;
  }
  return grad;
}

void CurvedFoldingBiasSignObjective::updateHessianIJV(const Eigen::VectorXd& x) {
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

		double t2 = ep_0_t-1.0;
		double t3 = ep_b_t-1.0;
		double t4 = t2*v2_z;
		double t5 = ep_f_t-1.0;
		double t6 = t2*v2_y;
		double t7 = ep_f_t*ep_f_v1_y;
		double t10 = ep_0_t*v1_y;
		double t23 = ep_f_v2_y*t5;
		double t8 = t6+t7-t10-t23;
		double t9 = ep_b_t*ep_b_v1_y;
		double t11 = t2*v2_x;
		double t12 = ep_b_t*ep_b_v1_x;
		double t19 = ep_0_t*v1_x;
		double t25 = ep_b_v2_x*t3;
		double t13 = t11+t12-t19-t25;
		double t14 = ep_f_t*ep_f_v1_z;
		double t17 = ep_0_t*v1_z;
		double t21 = ep_f_v2_z*t5;
		double t15 = t4+t14-t17-t21;
		double t16 = ep_b_t*ep_b_v1_z;
		double t18 = ep_f_t*ep_f_v1_x;
		double t27 = ep_b_v2_y*t3;
		double t20 = t6+t9-t10-t27;
		double t31 = ep_b_v2_z*t3;
		double t22 = t4+t16-t17-t31;
		double t24 = v1_z-v2_z;
		double t26 = t8*t13;
		double t32 = ep_f_v2_x*t5;
		double t28 = t11+t18-t19-t32;
		double t29 = v1_y-v2_y;
		double t30 = t13*t15;
		double t42 = t22*t28;
		double t33 = t30-t42;
		double t34 = t29*t33;
		double t35 = v1_x-v2_x;
		double t36 = t8*t22;
		double t40 = t15*t20;
		double t37 = t36-t40;
		double t38 = t35*t37;
		double t43 = t20*t28;
		double t44 = t26-t43;
		double t45 = t24*t44;
		double t46 = t34+t38-t45;
		double t47 = alpha*t46;
		double t39 = -sinh(t47);
		double t48 = t39*t39;
		double t41 = alpha*(t8*(t4+t16-ep_b_v2_z*t3-ep_0_t*v1_z)-t15*(t6+t9-ep_b_v2_y*t3-ep_0_t*v1_y))-alpha*t48*1.0/pow(cosh(alpha*(t34+t38-t24*(t26-t20*(t11+t18-ep_f_v2_x*t5-ep_0_t*v1_x)))),2.0)*(t36-t40);
		double t49 = cosh(t47);
		double t50 = 1.0/(t49*t49);
		double t51 = alpha*t37;
		double t55 = alpha*t48*t50*(t36-t40);
		double t52 = t51-t55;
		double t53 = alpha*t33;
		double t72 = alpha*t33*t48*t50;
		double t54 = t53-t72;
		double t56 = alpha*t44;
		double t73 = alpha*t44*t48*t50;
		double t57 = t56-t73;
		double t58 = w1_z-w2_z;
		double t59 = w1_y-w2_y;
		double t60 = t33*t59;
		double t61 = w1_x-w2_x;
		double t62 = t37*t61;
		double t65 = t44*t58;
		double t63 = t60+t62-t65;
		double t66 = alpha*t63;
		double t64 = -sinh(t66);
		double t67 = t64*t64;
		double t68 = cosh(t66);
		double t69 = 1.0/(t68*t68);
		double t76 = alpha*t33*t67*t69;
		double t70 = t53-t76;
		double t77 = alpha*t44*t67*t69;
		double t71 = t56-t77;
		double t74 = t54*t54;
		double t75 = t74*2.0;
		double t84 = alpha*t37*t67*t69;
		double t78 = t51-t84;
		double t79 = t54*t70*2.0;
		double t80 = t54*t57*2.0;
		double t81 = t57*(t51-t55)*2.0;
		double t82 = t57*t57;
		double t83 = t82*2.0;
		double t85 = t57*t71*2.0;
		double t86 = t51-t55;
		double t87 = t51-t55;
		double t88 = t54*(t51-t55)*2.0;
		double t89 = t71*(t51-t55)*2.0;
		double t90 = t70*(t51-t55)*2.0;
		double t91 = t54*t71*2.0;
		double t93 = alpha*t67*t69*(t36-t40);
		double t95 = t51-t93;
		double t92 = t54*t95*2.0;
		double t94 = t57*t70*2.0;
		double t96 = t57*(t51-t93)*2.0;
		double t97 = t51-t93;
		double t98 = t51-t93;
		double t99 = t70*(t51-t93)*2.0;
		double t100 = t70*t70;
		double t101 = t100*2.0;
		double t102 = t70*t71*2.0;
		double t103 = t71*(t51-t93)*2.0;
		double t104 = t71*t71;
		double t105 = t104*2.0;
		double t106 = (t51-t55)*(t51-t93)*2.0;
		double t107 = t54*(t51-t93)*2.0;
		double t108 = t51-t93;
		double t109 = t51-t93;

		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i, (t41*t41)*2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i+vnum, t54*(t51-alpha*t48*t50*(t36-t40))*2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v1_i+2*vnum, t57*(t51-alpha*t37*t48*t50)*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i, (t52*t52)*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i+vnum, t52*t54*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,v2_i+2*vnum, t81));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i, (t51-t55)*(t51-alpha*t67*1.0/pow(cosh(-alpha*t63),2.0)*(t36-t40))*2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i+vnum, t90));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w1_i+2*vnum, t52*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i, t52*t78*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i+vnum, t52*t70*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i,w2_i+2*vnum, t89));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i, t88));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i+vnum, t75));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v1_i+2*vnum, t54*t57*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i, t52*t54*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i+vnum, -t75));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,v2_i+2*vnum, t80));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i, t92));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i+vnum, t79));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w1_i+2*vnum, t54*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i, t54*t78*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i+vnum, -t79));
		IJV.push_back(Eigen::Triplet<double>(v1_i+vnum,w2_i+2*vnum, t91));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i, t52*t57*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i+vnum, -t80));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v1_i+2*vnum, t83));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i, t81));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i+vnum, t80));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,v2_i+2*vnum, -t83));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i, t57*t78*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i+vnum, t57*t70*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w1_i+2*vnum, t85));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i, t57*(t51-alpha*t67*t69*(t36-t40))*2.0));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i+vnum, t94));
		IJV.push_back(Eigen::Triplet<double>(v1_i+2*vnum,w2_i+2*vnum, -t85));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i, (t86*t86)*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i+vnum, t52*t54*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v1_i+2*vnum, t81));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i, (t87*t87)*2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i+vnum, t88));
		IJV.push_back(Eigen::Triplet<double>(v2_i,v2_i+2*vnum, t52*t57*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i, t52*t78*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i+vnum, t52*t70*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w1_i+2*vnum, t89));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i, (t51-alpha*t67*t69*(t36-t40))*(t51-t55)*2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i+vnum, t90));
		IJV.push_back(Eigen::Triplet<double>(v2_i,w2_i+2*vnum, t52*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i, t52*t54*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i+vnum, -t75));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v1_i+2*vnum, t80));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i, t88));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i+vnum, t75));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,v2_i+2*vnum, -t80));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i, t54*t78*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i+vnum, -t79));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w1_i+2*vnum, t91));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i, t92));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i+vnum, t79));
		IJV.push_back(Eigen::Triplet<double>(v2_i+vnum,w2_i+2*vnum, -t91));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i, t81));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i+vnum, t80));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v1_i+2*vnum, -t83));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i, t52*t57*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i+vnum, -t80));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,v2_i+2*vnum, t83));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i, t96));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i+vnum, t94));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w1_i+2*vnum, -t85));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i, t57*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i+vnum, -t94));
		IJV.push_back(Eigen::Triplet<double>(v2_i+2*vnum,w2_i+2*vnum, t85));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i, t106));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i+vnum, t107));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v1_i+2*vnum, t57*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i, t52*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i+vnum, -t92));
		IJV.push_back(Eigen::Triplet<double>(w1_i,v2_i+2*vnum, t96));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i, (t97*t97)*2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i+vnum, t99));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w1_i+2*vnum, t71*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i, (t98*t98)*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i+vnum, t70*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i,w2_i+2*vnum, t103));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i, t90));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i+vnum, t79));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v1_i+2*vnum, -t94));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i, t52*t70*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i+vnum, -t79));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,v2_i+2*vnum, t94));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i, t99));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i+vnum, t101));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w1_i+2*vnum, t70*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i, t70*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i+vnum, -t101));
		IJV.push_back(Eigen::Triplet<double>(w1_i+vnum,w2_i+2*vnum, t102));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i, t52*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i+vnum, -t91));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v1_i+2*vnum, t85));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i, t89));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i+vnum, t91));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,v2_i+2*vnum, -t85));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i, t71*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i+vnum, -t102));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w1_i+2*vnum, t105));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i, t103));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i+vnum, t102));
		IJV.push_back(Eigen::Triplet<double>(w1_i+2*vnum,w2_i+2*vnum, -t105));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i, t52*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i+vnum, -t92));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v1_i+2*vnum, t96));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i, t106));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i+vnum, t107));
		IJV.push_back(Eigen::Triplet<double>(w2_i,v2_i+2*vnum, t57*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i, (t108*t108)*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i+vnum, t70*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w1_i+2*vnum, t103));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i, (t109*t109)*2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i+vnum, t99));
		IJV.push_back(Eigen::Triplet<double>(w2_i,w2_i+2*vnum, t71*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i, t52*t70*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i+vnum, -t79));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v1_i+2*vnum, t94));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i, t90));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i+vnum, t79));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,v2_i+2*vnum, -t94));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i, t70*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i+vnum, -t101));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w1_i+2*vnum, t102));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i, t99));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i+vnum, t101));
		IJV.push_back(Eigen::Triplet<double>(w2_i+vnum,w2_i+2*vnum, -t102));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i, t89));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i+vnum, t91));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v1_i+2*vnum, -t85));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i, t52*t71*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i+vnum, -t91));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,v2_i+2*vnum, t85));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i, t103));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i+vnum, t102));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w1_i+2*vnum, -t105));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i, t71*t95*-2.0));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i+vnum, -t102));
		IJV.push_back(Eigen::Triplet<double>(w2_i+2*vnum,w2_i+2*vnum, t105));


  }
}