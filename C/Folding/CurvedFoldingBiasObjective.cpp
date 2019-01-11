#include "CurvedFoldingBiasObjective.h"

double CurvedFoldingBiasObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	int vnum = x.rows()/3;

	int h_cnt = 0;
	//for (int fold_i = 0; fold_i < curvedFoldBiases.size; fold_i++) {
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
		double t13 = ep_b_v2_x*t2;
		double t17 = ep_b_t*ep_b_v1_x;
		double t8 = t5+t7-t13-t17;
		double t9 = ep_f_v2_y*t6;
		double t11 = ep_b_v2_y*t2;
		double t26 = ep_f_t*ep_f_v1_y;
		double t10 = t3+t9-t11-t26;
		double t12 = ep_b_v2_y*t4;
		double t14 = ep_f_v2_x*t6;
		double t15 = ep_b_v2_z*t2;
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
		double t30 = (t18*(t15-t21+t22-ep_b_v2_z*t4)-t8*t20)*(-t3+t11+v2_y)+(t16-(t3+t12-ep_b_t*ep_b_v1_y-ep_b_v2_y*t2)*(t5+t14-ep_f_t*ep_f_v1_x-ep_b_v2_x*t2))*(t15+v1_z-ep_0_t*ep_0_v1_z)+(t8*t20-t18*(t15+t22-ep_0_t*ep_0_v1_z-ep_b_v2_z*t4))*(-t3+t11+v1_y)-(t16-t18*t24)*(t15+v2_z-ep_0_t*ep_0_v1_z)+t29*(-t5+t13+v1_x)-t29*(-t5+t13+v2_x);
		e += t30*t30;
  }
  return e;
}