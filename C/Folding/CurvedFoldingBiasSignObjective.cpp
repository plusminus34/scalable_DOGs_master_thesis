#include "CurvedFoldingBiasSignObjective.h"

void CurvedFoldingBiasSignObjective::reset_folds() {
	curvedFoldBiases.clear();
	//IJV.resize(0);
	//is_H_cached = false;	
}

void CurvedFoldingBiasSignObjective::add_fold_bias(const CurvedFoldBias& foldBias) {
	curvedFoldBiases.push_back(foldBias);
	//IJV.resize(curvedFoldBiases.size()*144);
	//is_H_cached = false; // recompute H caching next time
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
  // empty on purpose for now
}