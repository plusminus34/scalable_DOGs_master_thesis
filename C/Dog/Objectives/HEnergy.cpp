#include "HEnergy.h"

double HEnergy::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  int vnum = x.rows()/3;
    
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double t2 = p0_x-pxb_x;
		double t3 = p0_y-pxb_y;
		double t4 = p0_z-pxb_z;
		double t5 = p0_x-pxf_x;
		double t6 = p0_y-pxf_y;
		double t7 = p0_z-pxf_z;
		double t8 = t2*t2;
		double t9 = t3*t3;
		double t10 = t4*t4;
		double t11 = t8+t9+t10;
		double t12 = t5*t5;
		double t13 = t6*t6;
		double t14 = t7*t7;
		double t15 = t12+t13+t14;
		double t17 = 1.0/sqrt(t11);
		double t18 = 1.0/sqrt(t15);
		double t16 = t2*t17+t5*t18;
		double t19 = t3*t17+t6*t18;
		double t20 = t4*t17+t7*t18;
		double t21 = p0_x-pyb_x;
		double t22 = p0_y-pyb_y;
		double t23 = p0_z-pyb_z;
		double t24 = p0_x-pyf_x;
		double t25 = p0_y-pyf_y;
		double t26 = p0_z-pyf_z;
		double t27 = t21*t21;
		double t28 = t22*t22;
		double t29 = t23*t23;
		double t30 = t27+t28+t29;
		double t31 = t24*t24;
		double t32 = t25*t25;
		double t33 = t26*t26;
		double t34 = t31+t32+t33;
		double t36 = 1.0/sqrt(t30);
		double t37 = 1.0/sqrt(t34);
		double t35 = t21*t36+t24*t37;
		double t38 = t22*t36+t25*t37;
		double t39 = t23*t36+t26*t37;
		double t40 = sqrt(t16*t16+t19*t19+t20*t20)/(sqrt(t11)+sqrt(t15))+sqrt(t35*t35+t38*t38+t39*t39)/(sqrt(t30)+sqrt(t34));
		e +=  t40*t40;
  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));
    
	double t2 = p0_x-pxb_x;
	double t3 = p0_y-pxb_y;
	double t4 = p0_z-pxb_z;
	double t5 = p0_x-pxf_x;
	double t6 = p0_y-pxf_y;
	double t7 = p0_z-pxf_z;
	double t8 = t2*t2;
	double t9 = t3*t3;
	double t10 = t4*t4;
	double t11 = t8+t9+t10;
	double t12 = t5*t5;
	double t13 = t6*t6;
	double t14 = t7*t7;
	double t15 = t12+t13+t14;
	double t17 = 1.0/sqrt(t11);
	double t18 = 1.0/sqrt(t15);
	double t16 = t2*t17+t5*t18;
	double t19 = t3*t17+t6*t18;
	double t20 = t4*t17+t7*t18;
	e +=  (1.0/pow(sqrt(t11)+sqrt(t15),2.0)*(t16*t16+t19*t19+t20*t20));
  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  return e;
}

Eigen::VectorXd HEnergy::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;

  #pragma clang loop vectorize(enable)
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
        //local_grad.setZero();

        //int p_0_i = i*s+j, p_xf_i = i*s+j+1,p_xb_i = i*s+j-1,p_yf_i = (i+1)*s+j, p_yb_i = (i-1)*s+j;
        int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
        const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
        const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

        p0_z +=1e-17; // avoid dividing by zero

        /*
        Eigen::Vector3d e1b(p0_x-pxb_x,p0_y-pxb_y,p0_z-pxb_z);
		Eigen::Vector3d e1(p0_x-pxf_x,p0_y-pxf_y,p0_z-pxf_z);
		Eigen::Vector3d e2b(p0_x-pyb_x,p0_y-pyb_y,p0_z-pyb_z);
		Eigen::Vector3d e2(p0_x-pyf_x,p0_y-pyf_y,p0_z-pyf_z);
		*/

		double t2 = p0_x-pxb_x;
		double t3 = p0_y-pxb_y;
		double t4 = p0_z-pxb_z;
		double t5 = p0_x-pxf_x;
		double t6 = p0_y-pxf_y;
		double t7 = p0_z-pxf_z;
		double t8 = t2*t2;
		double t9 = t3*t3;
		double t10 = t4*t4;
		double t11 = t8+t9+t10;
		double t12 = t5*t5;
		double t13 = t6*t6;
		double t14 = t7*t7;
		double t15 = t12+t13+t14;
		double t17 = 1.0/sqrt(t11);
		double t18 = 1.0/sqrt(t15);
		double t44 = t2*t17;
		double t45 = t5*t18;
		double t16 = t44+t45;
		double t47 = t3*t17;
		double t48 = t6*t18;
		double t19 = t47+t48;
		double t50 = t4*t17;
		double t51 = t7*t18;
		double t20 = t50+t51;
		double t21 = p0_x-pyb_x;
		double t22 = p0_y-pyb_y;
		double t23 = p0_z-pyb_z;
		double t24 = p0_x-pyf_x;
		double t25 = p0_y-pyf_y;
		double t26 = p0_z-pyf_z;
		double t27 = t21*t21;
		double t28 = t22*t22;
		double t29 = t23*t23;
		double t30 = t27+t28+t29;
		double t31 = t24*t24;
		double t32 = t25*t25;
		double t33 = t26*t26;
		double t34 = t31+t32+t33;
		double t36 = 1.0/sqrt(t30);
		double t37 = 1.0/sqrt(t34);
		double t58 = t21*t36;
		double t59 = t24*t37;
		double t35 = t58+t59;
		double t61 = t22*t36;
		double t62 = t25*t37;
		double t38 = t61+t62;
		double t64 = t23*t36;
		double t65 = t26*t37;
		double t39 = t64+t65;
		double t40 = p0_x*2.0;
		double t41 = sqrt(t11);
		double t42 = sqrt(t15);
		double t43 = t41+t42;
		double t46 = t16*t16;
		double t49 = t19*t19;
		double t52 = t20*t20;
		double t53 = t46+t49+t52;
		double t54 = sqrt(t53);
		double t55 = sqrt(t30);
		double t56 = sqrt(t34);
		double t57 = t55+t56;
		double t60 = t35*t35;
		double t63 = t38*t38;
		double t66 = t39*t39;
		double t67 = t60+t63+t66;
		double t68 = sqrt(t67);
		double t69 = 1.0/t43;
		double t70 = pxb_x*2.0;
		double t71 = t40-t70;
		double t72 = pxf_x*2.0;
		double t73 = t40-t72;
		double t74 = 1.0/pow(t11,3.0/2.0);
		double t75 = 1.0/pow(t15,3.0/2.0);
		double t76 = 1.0/t57;
		double t77 = pyb_x*2.0;
		double t78 = t40-t77;
		double t79 = pyf_x*2.0;
		double t80 = t40-t79;
		double t81 = 1.0/pow(t30,3.0/2.0);
		double t82 = 1.0/pow(t34,3.0/2.0);
		double t83 = t54*t69;
		double t84 = t68*t76;
		double t85 = t83+t84;
		double t86 = p0_y*2.0;
		double t87 = 1.0/(t43*t43);
		double t88 = 1.0/(t57*t57);
		double t89 = pxb_y*2.0;
		double t90 = t86-t89;
		double t91 = pxf_y*2.0;
		double t92 = t86-t91;
		double t93 = 1.0/sqrt(t53); // this is infinity for straight-x
		double t94 = pyb_y*2.0;
		double t95 = t86-t94;
		double t96 = pyf_y*2.0;
		double t97 = t86-t96;
		double t98 = 1.0/sqrt(t67); // this is infinity for straight-y
		double t99 = p0_z*2.0;
		double t100 = pxb_z*2.0;
		double t101 = t99-t100;
		double t102 = pxf_z*2.0;
		double t103 = t99-t102;
		double t104 = pyb_z*2.0;
		double t105 = t99-t104;
		double t106 = pyf_z*2.0;
		double t107 = t99-t106;

		/*
		bool x_straight = false, y_straight = false;
		double eps = 1e-12;
		if (e1.cross(e1b).squaredNorm() < eps) {
			//cout << "straight bnd!" << endl;
			x_straight = true;
		}
		if (e2.cross(e2b).squaredNorm() < eps) {
			//cout << "straight bnd!" << endl;
			y_straight = true;
		}
		//cout << "before grad(p_0_i) = " << grad(p_0_i) << endl;
		//cout << "before grad(p_0_i+v_num) = " << grad(p_0_i+v_num) << endl;
		//cout << "before grad(p_0_i+2*+v_num) = " << grad(p_0_i+2*v_num) << endl;

		//grad(p_0_i) +=  (t85*(t54*t87*(t17*t71*(1.0/2.0)+t18*t73*(1.0/2.0))+t68*t88*(t36*t78*(1.0/2.0)+t37*t80*(1.0/2.0))+t69*t93*(t16*(t17+t18-t2*t71*t74*(1.0/2.0)-t5*t73*t75*(1.0/2.0))*-2.0+t19*(t3*t71*t74*(1.0/2.0)+t6*t73*t75*(1.0/2.0))*2.0+t20*(t4*t71*t74*(1.0/2.0)+t7*t73*t75*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t35*(t36+t37-t21*t78*t81*(1.0/2.0)-t24*t80*t82*(1.0/2.0))*-2.0+t38*(t22*t78*t81*(1.0/2.0)+t25*t80*t82*(1.0/2.0))*2.0+t39*(t23*t78*t81*(1.0/2.0)+t26*t80*t82*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
		//grad(p_0_i+v_num) +=  (t85*(t54*t87*(t17*t90*(1.0/2.0)+t18*t92*(1.0/2.0))+t68*t88*(t36*t95*(1.0/2.0)+t37*t97*(1.0/2.0))+t69*t93*(t19*(t17+t18-t3*t74*t90*(1.0/2.0)-t6*t75*t92*(1.0/2.0))*-2.0+t16*(t2*t74*t90*(1.0/2.0)+t5*t75*t92*(1.0/2.0))*2.0+t20*(t4*t74*t90*(1.0/2.0)+t7*t75*t92*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t38*(t36+t37-t22*t81*t95*(1.0/2.0)-t25*t82*t97*(1.0/2.0))*-2.0+t35*(t21*t81*t95*(1.0/2.0)+t24*t82*t97*(1.0/2.0))*2.0+t39*(t23*t81*t95*(1.0/2.0)+t26*t82*t97*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
		//grad(p_0_i+2*v_num) +=  (t85*(t54*t87*(t17*t101*(1.0/2.0)+t18*t103*(1.0/2.0))+t68*t88*(t36*t105*(1.0/2.0)+t37*t107*(1.0/2.0))+t69*t93*(t20*(t17+t18-t4*t74*t101*(1.0/2.0)-t7*t75*t103*(1.0/2.0))*-2.0+t16*(t2*t74*t101*(1.0/2.0)+t5*t75*t103*(1.0/2.0))*2.0+t19*(t3*t74*t101*(1.0/2.0)+t6*t75*t103*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t39*(t36+t37-t23*t81*t105*(1.0/2.0)-t26*t82*t107*(1.0/2.0))*-2.0+t35*(t21*t81*t105*(1.0/2.0)+t24*t82*t107*(1.0/2.0))*2.0+t38*(t22*t81*t105*(1.0/2.0)+t25*t82*t107*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);

		//cout << "after grad(p_0_i) = " << grad(p_0_i) << endl;
		//cout << "after grad(p_0_i+v_num) = " << grad(p_0_i+v_num) << endl;
		//cout << "after grad(p_0_i+2*+v_num) = " << grad(p_0_i+2*v_num) << endl;

		cout << "1 line is " << t54*t87*(t17*t71*(1.0/2.0)+t18*t73*(1.0/2.0)) << endl;
			cout << "2 line is " << t68*t88*(t36*t78*(1.0/2.0)+t37*t80*(1.0/2.0)) << endl;
			cout << "3 lines is " << t69*t93*(t16*(t17+t18-t2*t71*t74*(1.0/2.0)-t5*t73*t75*(1.0/2.0))*-2.0+t19*(t3*t71*t74*(1.0/2.0)+t6*t73*t75*(1.0/2.0))*2.0+t20*(t4*t71*t74*(1.0/2.0)+t7*t73*t75*(1.0/2.0))*2.0)*(1.0/2.0) << endl;
			cout << "4 lines is " << t76*t98*(t35*(t36+t37-t21*t78*t81*(1.0/2.0)-t24*t80*t82*(1.0/2.0))*-2.0+t38*(t22*t78*t81*(1.0/2.0)+t25*t80*t82*(1.0/2.0))*2.0+t39*(t23*t78*t81*(1.0/2.0)+t26*t80*t82*(1.0/2.0))*2.0)*(1.0/2.0) << endl;

		if (!x_straight) {
			grad(p_0_i) +=  (t85*
										(t54*t87*(t17*t71*(1.0/2.0)+t18*t73*(1.0/2.0))
										+t68*t88*(t36*t78*(1.0/2.0)+t37*t80*(1.0/2.0))
										+t69*t93*(t16*(t17+t18-t2*t71*t74*(1.0/2.0)-t5*t73*t75*(1.0/2.0))*-2.0+t19*(t3*t71*t74*(1.0/2.0)+t6*t73*t75*(1.0/2.0))*2.0+t20*(t4*t71*t74*(1.0/2.0)+t7*t73*t75*(1.0/2.0))*2.0)*(1.0/2.0)
										+t76*t98*(t35*(t36+t37-t21*t78*t81*(1.0/2.0)-t24*t80*t82*(1.0/2.0))*-2.0+t38*(t22*t78*t81*(1.0/2.0)+t25*t80*t82*(1.0/2.0))*2.0+t39*(t23*t78*t81*(1.0/2.0)+t26*t80*t82*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
			cout << "1 line is " << t54*t87*(t17*t71*(1.0/2.0)+t18*t73*(1.0/2.0)) << endl;
			cout << "2 line is " << t68*t88*(t36*t78*(1.0/2.0)+t37*t80*(1.0/2.0)) << endl;
			cout << "3 lines is " << t69*t93*(t16*(t17+t18-t2*t71*t74*(1.0/2.0)-t5*t73*t75*(1.0/2.0))*-2.0+t19*(t3*t71*t74*(1.0/2.0)+t6*t73*t75*(1.0/2.0))*2.0+t20*(t4*t71*t74*(1.0/2.0)+t7*t73*t75*(1.0/2.0))*2.0)*(1.0/2.0) << endl;
			cout << "4 lines is " << t76*t98*(t35*(t36+t37-t21*t78*t81*(1.0/2.0)-t24*t80*t82*(1.0/2.0))*-2.0+t38*(t22*t78*t81*(1.0/2.0)+t25*t80*t82*(1.0/2.0))*2.0+t39*(t23*t78*t81*(1.0/2.0)+t26*t80*t82*(1.0/2.0))*2.0)*(1.0/2.0) << endl;

			grad(p_xb_i+0) +=  (t85*(t69*t93*(t16*(t17-t2*t71*t74*(1.0/2.0))*-2.0+t3*t19*t71*t74+t4*t20*t71*t74)*(1.0/2.0)+t17*t54*t71*t87*(1.0/2.0))*2.0);
			grad(p_xb_i+v_num) +=  (t85*(t69*t93*(t19*(t17-t3*t74*t90*(1.0/2.0))*-2.0+t2*t16*t74*t90+t4*t20*t74*t90)*(1.0/2.0)+t17*t54*t87*t90*(1.0/2.0))*2.0);
			grad(p_xb_i+2*v_num) +=  (t85*(t69*t93*(t20*(t17-t4*t74*t101*(1.0/2.0))*-2.0+t2*t16*t74*t101+t3*t19*t74*t101)*(1.0/2.0)+t17*t54*t87*t101*(1.0/2.0))*2.0);
			grad(p_xf_i) +=  (t85*(t69*t93*(t16*(t18-t5*t73*t75*(1.0/2.0))*-2.0+t6*t19*t73*t75+t7*t20*t73*t75)*(1.0/2.0)+t18*t54*t73*t87*(1.0/2.0))*2.0);
			grad(p_xf_i+v_num) +=  (t85*(t69*t93*(t19*(t18-t6*t75*t92*(1.0/2.0))*-2.0+t5*t16*t75*t92+t7*t20*t75*t92)*(1.0/2.0)+t18*t54*t87*t92*(1.0/2.0))*2.0);
			grad(p_xf_i+2*v_num) +=  (t85*(t69*t93*(t20*(t18-t7*t75*t103*(1.0/2.0))*-2.0+t5*t16*t75*t103+t6*t19*t75*t103)*(1.0/2.0)+t18*t54*t87*t103*(1.0/2.0))*2.0);
		}
		if (!y_straight) {
			grad(p_yb_i) +=  (t85*(t76*t98*(t35*(t36-t21*t78*t81*(1.0/2.0))*-2.0+t22*t38*t78*t81+t23*t39*t78*t81)*(1.0/2.0)+t36*t68*t78*t88*(1.0/2.0))*2.0);
			grad(p_yb_i+v_num) +=  (t85*(t76*t98*(t38*(t36-t22*t81*t95*(1.0/2.0))*-2.0+t21*t35*t81*t95+t23*t39*t81*t95)*(1.0/2.0)+t36*t68*t88*t95*(1.0/2.0))*2.0);
			grad(p_yb_i+2*v_num) +=  (t85*(t76*t98*(t39*(t36-t23*t81*t105*(1.0/2.0))*-2.0+t21*t35*t81*t105+t22*t38*t81*t105)*(1.0/2.0)+t36*t68*t88*t105*(1.0/2.0))*2.0);
			grad(p_yf_i) +=  (t85*(t76*t98*(t35*(t37-t24*t80*t82*(1.0/2.0))*-2.0+t25*t38*t80*t82+t26*t39*t80*t82)*(1.0/2.0)+t37*t68*t80*t88*(1.0/2.0))*2.0);
			grad(p_yf_i+v_num) +=  (t85*(t76*t98*(t38*(t37-t25*t82*t97*(1.0/2.0))*-2.0+t24*t35*t82*t97+t26*t39*t82*t97)*(1.0/2.0)+t37*t68*t88*t97*(1.0/2.0))*2.0);
			grad(p_yf_i+2*v_num) +=  (t85*(t76*t98*(t39*(t37-t26*t82*t107*(1.0/2.0))*-2.0+t24*t35*t82*t107+t25*t38*t82*t107)*(1.0/2.0)+t37*t68*t88*t107*(1.0/2.0))*2.0);
		}*/
		grad(p_0_i) +=  (t85*(t54*t87*(t17*t71*(1.0/2.0)+t18*t73*(1.0/2.0))+t68*t88*(t36*t78*(1.0/2.0)+t37*t80*(1.0/2.0))+t69*t93*(t16*(t17+t18-t2*t71*t74*(1.0/2.0)-t5*t73*t75*(1.0/2.0))*-2.0+t19*(t3*t71*t74*(1.0/2.0)+t6*t73*t75*(1.0/2.0))*2.0+t20*(t4*t71*t74*(1.0/2.0)+t7*t73*t75*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t35*(t36+t37-t21*t78*t81*(1.0/2.0)-t24*t80*t82*(1.0/2.0))*-2.0+t38*(t22*t78*t81*(1.0/2.0)+t25*t80*t82*(1.0/2.0))*2.0+t39*(t23*t78*t81*(1.0/2.0)+t26*t80*t82*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
		grad(p_0_i+v_num) +=  (t85*(t54*t87*(t17*t90*(1.0/2.0)+t18*t92*(1.0/2.0))+t68*t88*(t36*t95*(1.0/2.0)+t37*t97*(1.0/2.0))+t69*t93*(t19*(t17+t18-t3*t74*t90*(1.0/2.0)-t6*t75*t92*(1.0/2.0))*-2.0+t16*(t2*t74*t90*(1.0/2.0)+t5*t75*t92*(1.0/2.0))*2.0+t20*(t4*t74*t90*(1.0/2.0)+t7*t75*t92*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t38*(t36+t37-t22*t81*t95*(1.0/2.0)-t25*t82*t97*(1.0/2.0))*-2.0+t35*(t21*t81*t95*(1.0/2.0)+t24*t82*t97*(1.0/2.0))*2.0+t39*(t23*t81*t95*(1.0/2.0)+t26*t82*t97*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
		grad(p_0_i+2*v_num) +=  (t85*(t54*t87*(t17*t101*(1.0/2.0)+t18*t103*(1.0/2.0))+t68*t88*(t36*t105*(1.0/2.0)+t37*t107*(1.0/2.0))+t69*t93*(t20*(t17+t18-t4*t74*t101*(1.0/2.0)-t7*t75*t103*(1.0/2.0))*-2.0+t16*(t2*t74*t101*(1.0/2.0)+t5*t75*t103*(1.0/2.0))*2.0+t19*(t3*t74*t101*(1.0/2.0)+t6*t75*t103*(1.0/2.0))*2.0)*(1.0/2.0)+t76*t98*(t39*(t36+t37-t23*t81*t105*(1.0/2.0)-t26*t82*t107*(1.0/2.0))*-2.0+t35*(t21*t81*t105*(1.0/2.0)+t24*t82*t107*(1.0/2.0))*2.0+t38*(t22*t81*t105*(1.0/2.0)+t25*t82*t107*(1.0/2.0))*2.0)*(1.0/2.0))*-2.0);
		grad(p_xb_i+0) +=  (t85*(t69*t93*(t16*(t17-t2*t71*t74*(1.0/2.0))*-2.0+t3*t19*t71*t74+t4*t20*t71*t74)*(1.0/2.0)+t17*t54*t71*t87*(1.0/2.0))*2.0);
		grad(p_xb_i+v_num) +=  (t85*(t69*t93*(t19*(t17-t3*t74*t90*(1.0/2.0))*-2.0+t2*t16*t74*t90+t4*t20*t74*t90)*(1.0/2.0)+t17*t54*t87*t90*(1.0/2.0))*2.0);
		grad(p_xb_i+2*v_num) +=  (t85*(t69*t93*(t20*(t17-t4*t74*t101*(1.0/2.0))*-2.0+t2*t16*t74*t101+t3*t19*t74*t101)*(1.0/2.0)+t17*t54*t87*t101*(1.0/2.0))*2.0);
		grad(p_xf_i) +=  (t85*(t69*t93*(t16*(t18-t5*t73*t75*(1.0/2.0))*-2.0+t6*t19*t73*t75+t7*t20*t73*t75)*(1.0/2.0)+t18*t54*t73*t87*(1.0/2.0))*2.0);
		grad(p_xf_i+v_num) +=  (t85*(t69*t93*(t19*(t18-t6*t75*t92*(1.0/2.0))*-2.0+t5*t16*t75*t92+t7*t20*t75*t92)*(1.0/2.0)+t18*t54*t87*t92*(1.0/2.0))*2.0);
		grad(p_xf_i+2*v_num) +=  (t85*(t69*t93*(t20*(t18-t7*t75*t103*(1.0/2.0))*-2.0+t5*t16*t75*t103+t6*t19*t75*t103)*(1.0/2.0)+t18*t54*t87*t103*(1.0/2.0))*2.0);
		grad(p_yb_i) +=  (t85*(t76*t98*(t35*(t36-t21*t78*t81*(1.0/2.0))*-2.0+t22*t38*t78*t81+t23*t39*t78*t81)*(1.0/2.0)+t36*t68*t78*t88*(1.0/2.0))*2.0);
		grad(p_yb_i+v_num) +=  (t85*(t76*t98*(t38*(t36-t22*t81*t95*(1.0/2.0))*-2.0+t21*t35*t81*t95+t23*t39*t81*t95)*(1.0/2.0)+t36*t68*t88*t95*(1.0/2.0))*2.0);
		grad(p_yb_i+2*v_num) +=  (t85*(t76*t98*(t39*(t36-t23*t81*t105*(1.0/2.0))*-2.0+t21*t35*t81*t105+t22*t38*t81*t105)*(1.0/2.0)+t36*t68*t88*t105*(1.0/2.0))*2.0);
		grad(p_yf_i) +=  (t85*(t76*t98*(t35*(t37-t24*t80*t82*(1.0/2.0))*-2.0+t25*t38*t80*t82+t26*t39*t80*t82)*(1.0/2.0)+t37*t68*t80*t88*(1.0/2.0))*2.0);
		grad(p_yf_i+v_num) +=  (t85*(t76*t98*(t38*(t37-t25*t82*t97*(1.0/2.0))*-2.0+t24*t35*t82*t97+t26*t39*t82*t97)*(1.0/2.0)+t37*t68*t88*t97*(1.0/2.0))*2.0);
		grad(p_yf_i+2*v_num) +=  (t85*(t76*t98*(t39*(t37-t26*t82*t107*(1.0/2.0))*-2.0+t24*t35*t82*t107+t25*t38*t82*t107)*(1.0/2.0)+t37*t68*t88*t107*(1.0/2.0))*2.0);
    //}
  }


  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	double t2 = p0_x-pxb_x;
	double t3 = p0_y-pxb_y;
	double t4 = p0_z-pxb_z;
	double t5 = p0_x-pxf_x;
	double t6 = p0_y-pxf_y;
	double t7 = p0_z-pxf_z;
	double t8 = t2*t2;
	double t9 = t3*t3;
	double t10 = t4*t4;
	double t11 = t8+t9+t10;
	double t12 = p0_x*2.0;
	double t13 = t5*t5;
	double t14 = t6*t6;
	double t15 = t7*t7;
	double t16 = t13+t14+t15;
	double t17 = pxb_x*2.0;
	double t18 = t12-t17;
	double t19 = 1.0/pow(t11,3.0/2.0);
	double t20 = pxf_x*2.0;
	double t21 = t12-t20;
	double t22 = 1.0/pow(t16,3.0/2.0);
	double t23 = 1.0/sqrt(t11);
	double t24 = 1.0/sqrt(t16);
	double t25 = sqrt(t11);
	double t26 = sqrt(t16);
	double t27 = t25+t26;
	double t28 = t2*t23;
	double t29 = t5*t24;
	double t30 = t28+t29;
	double t31 = t3*t23;
	double t32 = t6*t24;
	double t33 = t31+t32;
	double t34 = t4*t23;
	double t35 = t7*t24;
	double t36 = t34+t35;
	double t37 = 1.0/(t27*t27);
	double t38 = p0_y*2.0;
	double t39 = pxb_y*2.0;
	double t40 = t38-t39;
	double t41 = pxf_y*2.0;
	double t42 = t38-t41;
	double t43 = 1.0/(t27*t27*t27);
	double t44 = t30*t30;
	double t45 = t33*t33;
	double t46 = t36*t36;
	double t47 = t44+t45+t46;
	double t48 = p0_z*2.0;
	double t49 = pxb_z*2.0;
	double t50 = t48-t49;
	double t51 = pxf_z*2.0;
	double t52 = t48-t51;
	
	grad(p_0_i) += (-t37*(t30*(t23+t24-t2*t18*t19*(1.0/2.0)-t5*t21*t22*(1.0/2.0))*-2.0+t33*(t3*t18*t19*(1.0/2.0)+t6*t21*t22*(1.0/2.0))*2.0+t36*(t4*t18*t19*(1.0/2.0)+t7*t21*t22*(1.0/2.0))*2.0)-t43*t47*(t18*t23*(1.0/2.0)+t21*t24*(1.0/2.0))*2.0);
	grad(p_0_i+v_num) += (-t37*(t33*(t23+t24-t3*t19*t40*(1.0/2.0)-t6*t22*t42*(1.0/2.0))*-2.0+t30*(t2*t19*t40*(1.0/2.0)+t5*t22*t42*(1.0/2.0))*2.0+t36*(t4*t19*t40*(1.0/2.0)+t7*t22*t42*(1.0/2.0))*2.0)-t43*t47*(t23*t40*(1.0/2.0)+t24*t42*(1.0/2.0))*2.0);
	grad(p_0_i+2*v_num) += (-t37*(t36*(t23+t24-t4*t19*t50*(1.0/2.0)-t7*t22*t52*(1.0/2.0))*-2.0+t30*(t2*t19*t50*(1.0/2.0)+t5*t22*t52*(1.0/2.0))*2.0+t33*(t3*t19*t50*(1.0/2.0)+t6*t22*t52*(1.0/2.0))*2.0)-t43*t47*(t23*t50*(1.0/2.0)+t24*t52*(1.0/2.0))*2.0);
	grad(p_xb_i+0) += (t37*(t30*(t23-t2*t18*t19*(1.0/2.0))*-2.0+t3*t18*t19*t33+t4*t18*t19*t36)+t18*t23*t43*t47);
	grad(p_xb_i+v_num) += (t37*(t33*(t23-t3*t19*t40*(1.0/2.0))*-2.0+t2*t19*t30*t40+t4*t19*t36*t40)+t23*t40*t43*t47);
	grad(p_xb_i+2*v_num) += (t37*(t36*(t23-t4*t19*t50*(1.0/2.0))*-2.0+t2*t19*t30*t50+t3*t19*t33*t50)+t23*t43*t47*t50);
	grad(p_xf_i) += (t37*(t30*(t24-t5*t21*t22*(1.0/2.0))*-2.0+t6*t21*t22*t33+t7*t21*t22*t36)+t21*t24*t43*t47);
	grad(p_xf_i+v_num) += (t37*(t33*(t24-t6*t22*t42*(1.0/2.0))*-2.0+t5*t22*t30*t42+t7*t22*t36*t42)+t24*t42*t43*t47);
	grad(p_xf_i+2*v_num) += (t37*(t36*(t24-t7*t22*t52*(1.0/2.0))*-2.0+t5*t22*t30*t52+t6*t22*t33*t52)+t24*t43*t47*t52);

  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  //cout << "new grad.norm() = " << grad.norm() << endl; //exit(1);
  return grad;
}