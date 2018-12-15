#include "EqualDiagObjective.h"

double EqualDiagObjective::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  int vnum = x.rows()/3;
    
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double t2 = pxb_x-pyb_x;
		double t3 = pxb_y-pyb_y;
		double t4 = pxb_z-pyb_z;
		double t5 = pxb_x-pyf_x;
		double t6 = pxb_y-pyf_y;
		double t7 = pxb_z-pyf_z;
		double t9 = t2*t2;
		double t10 = t3*t3;
		double t11 = t4*t4;
		double t16 = t5*t5;
		double t17 = t6*t6;
		double t18 = t7*t7;
		double t8 = t9+t10+t11-t16-t17-t18;
		double t12 = pxf_x-pyb_x;
		double t13 = pxf_y-pyb_y;
		double t14 = pxf_z-pyb_z;
		double t23 = t12*t12;
		double t24 = t13*t13;
		double t25 = t14*t14;
		double t15 = t9+t10+t11-t23-t24-t25;
		double t19 = pxf_x-pyf_x;
		double t20 = pxf_y-pyf_y;
		double t21 = pxf_z-pyf_z;
		double t26 = t19*t19;
		double t27 = t20*t20;
		double t28 = t21*t21;
		double t22 = t16+t17+t18-t26-t27-t28;
		double t29 = t23+t24+t25-t26-t27-t28;
		e += t8*t8+t15*t15+t22*t22+t29*t29;
  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));
    
  	double t2 = pxb_x-pyf_x;
  	double t3 = pxb_y-pyf_y;
  	double t4 = pxb_z-pyf_z;
  	double t5 = pxf_x-pyf_x;
    double t6 = pxf_y-pyf_y;
    double t7 = pxf_z-pyf_z;
    double t8 = t2*t2+t3*t3+t4*t4-t5*t5-t6*t6-t7*t7;
	e += t8*t8;
  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  return e;
}

Eigen::VectorXd EqualDiagObjective::grad(const Eigen::VectorXd& x) const {
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
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  const double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

        double t2 = pxb_x-pyb_x;
		double t3 = pxb_y-pyb_y;
		double t4 = pxb_z-pyb_z;
		double t5 = pxf_x-pyb_x;
		double t6 = pxf_y-pyb_y;
		double t7 = pxf_z-pyb_z;
		double t8 = pxb_x*2.0;
		double t9 = pxb_x-pyf_x;
		double t10 = pxb_y-pyf_y;
		double t11 = pxb_z-pyf_z;
		double t12 = pxf_x-pyf_x;
		double t13 = pxf_y-pyf_y;
		double t14 = pxf_z-pyf_z;
		double t15 = pyb_x*2.0;
		double t16 = pyf_x*2.0;
		double t17 = t2*t2;
		double t18 = t3*t3;
		double t19 = t4*t4;
		double t20 = t9*t9;
		double t21 = t10*t10;
		double t22 = t11*t11;
		double t23 = t5*t5;
		double t24 = t6*t6;
		double t25 = t7*t7;
		double t26 = t17+t18+t19-t23-t24-t25;
		double t27 = pxb_y*2.0;
		double t28 = t12*t12;
		double t29 = t13*t13;
		double t30 = t14*t14;
		double t31 = t20+t21+t22-t28-t29-t30;
		double t32 = pyb_y*2.0;
		double t33 = pyf_y*2.0;
		double t34 = t17+t18+t19-t20-t21-t22;
		double t35 = pxb_z*2.0;
		double t36 = pyb_z*2.0;
		double t37 = pyf_z*2.0;
		double t38 = pxf_x*2.0;
		double t39 = t15-t16;
		double t40 = pxf_y*2.0;
		double t41 = t32-t33;
		double t42 = t23+t24+t25-t28-t29-t30;
		double t43 = pxf_z*2.0;
		double t44 = t36-t37;
		double t45 = t8-t15;
		double t46 = t27-t32;
		double t47 = t35-t36;
		double t48 = t8-t38;
		double t49 = t8-t16;
		double t50 = t16-t38;
		double t51 = t27-t40;
		double t52 = t27-t33;
		double t53 = t33-t40;
		double t54 = t35-t43;
		double t55 = t35-t37;
		double t56 = t37-t43;

		grad(p_xb_i+0) += t26*t45*2.0-t34*t39*2.0+t31*t49*2.0;
		grad(p_xb_i+v_num) += t26*t46*2.0-t34*t41*2.0+t31*t52*2.0;
		grad(p_xb_i+2*v_num) += t26*t47*2.0-t34*t44*2.0+t31*t55*2.0;
		grad(p_xf_i) += t31*t50*2.0-t39*t42*2.0+t26*(t15-t38)*2.0;
		grad(p_xf_i+v_num) += t41*t42*-2.0+t31*t53*2.0+t26*(t32-t40)*2.0;
		grad(p_xf_i+2*v_num) += t42*t44*-2.0+t31*t56*2.0+t26*(t36-t43)*2.0;
		grad(p_yb_i) += t26*t48*-2.0-t34*t45*2.0+t42*(t15-t38)*2.0;
		grad(p_yb_i+v_num) += t26*t51*-2.0-t34*t46*2.0+t42*(t32-t40)*2.0;
		grad(p_yb_i+2*v_num) += t26*t54*-2.0-t34*t47*2.0+t42*(t36-t43)*2.0;
		grad(p_yf_i) += t31*t48*-2.0+t34*t49*2.0-t42*t50*2.0;
		grad(p_yf_i+v_num) += t31*t51*-2.0+t34*t52*2.0-t42*t53*2.0;
		grad(p_yf_i+2*v_num) += t31*t54*-2.0+t34*t55*2.0-t42*t56*2.0;

  }


  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	double t2 = pxb_x-pyf_x;
	double t3 = pxb_y-pyf_y;
	double t4 = pxb_z-pyf_z;
	double t5 = pxf_x-pyf_x;
	double t6 = pxf_y-pyf_y;
	double t7 = pxf_z-pyf_z;
	double t8 = t2*t2;
	double t9 = t3*t3;
	double t10 = t4*t4;
	double t11 = t5*t5;
	double t12 = t6*t6;
	double t13 = t7*t7;
	double t14 = t8+t9+t10-t11-t12-t13;
	double t15 = pxb_x*2.0;
	double t16 = pxf_x*2.0;
	double t17 = pxb_y*2.0;
	double t18 = pxf_y*2.0;
	double t19 = pxb_z*2.0;
	double t20 = pxf_z*2.0;

	grad(p_xb_i+0) += t14*(pyf_x*2.0-t15)*-2.0;
	grad(p_xb_i+v_num) += t14*(pyf_y*2.0-t17)*-2.0;
	grad(p_xb_i+2*v_num) += t14*(pyf_z*2.0-t19)*-2.0;
	grad(p_xf_i) += t14*(pyf_x*2.0-t16)*2.0;
	grad(p_xf_i+v_num) += t14*(pyf_y*2.0-t18)*2.0;
	grad(p_xf_i+2*v_num) += t14*(pyf_z*2.0-t20)*2.0;
	grad(p_yf_i) +=  t14*(t15-t16)*-2.0;
	grad(p_yf_i+v_num) += t14*(t17-t18)*-2.0;
	grad(p_yf_i+2*v_num) += t14*(t19-t20)*-2.0;
	
  }
  //cout << "new grad.norm() = " << grad.norm() << endl; //exit(1);
  return grad;
}