#include "IsometryObjective.h"

IsometryObjective::IsometryObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {
	refL.resize(quadTop.E.rows()); refL.setZero();
}

void IsometryObjective::set_ref(const Eigen::VectorXd& x) {
	int vnum = x.rows()/3;

  	int h_cnt = 0;
  	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

		double t2 = p0_x-pxf_x;
		double t3 = p0_y-pxf_y;
		double t4 = p0_z-pxf_z;
		refL[h_cnt] = t2*t2+t3*t3+t4*t4; // squared length
		h_cnt++;
  }
}

double IsometryObjective::obj(const Eigen::VectorXd& x) {
	double e = 0;
	int vnum = x.rows()/3;

	int h_cnt = 0;
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

  		double t2 = p0_x-pxf_x;
  		double t3 = p0_y-pxf_y;
  		double t4 = p0_z-pxf_z;
  		double t5 = -l0+t2*t2+t3*t3+t4*t4;
  		e += t5*t5;

		h_cnt++;
  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  return e;
}

Eigen::VectorXd IsometryObjective::grad(const Eigen::VectorXd& x) {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0;

  #pragma clang loop vectorize(enable)
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

		double t2 = p0_x-pxf_x;
		double t3 = p0_y-pxf_y;
		double t4 = p0_z-pxf_z;
		double t5 = t2*t2;
		double t6 = t3*t3;
		double t7 = t4*t4;
		double t8 = -l0+t5+t6+t7;
		double t9 = p0_x*2.0;
		double t10 = pxf_x*2.0;
		double t11 = t9-t10;
		double t12 = t8*t11*2.0;
		double t13 = p0_y*2.0;
		double t14 = pxf_y*2.0;
		double t15 = t13-t14;
		double t16 = t8*t15*2.0;
		double t17 = p0_z*2.0;
		double t18 = pxf_z*2.0;
		double t19 = t17-t18;
		double t20 = t8*t19*2.0;

		grad(p_0_i) += t12;
        grad(p_0_i+v_num) += t16;
        grad(p_0_i+2*v_num) += t20;
        grad(p_xf_i) += -t12;
        grad(p_xf_i+v_num) += -t16;
        grad(p_xf_i+2*v_num) += -t20;

        h_cnt++;
  }
  return grad;
}

Eigen::SparseMatrix<double> IsometryObjective::hessian(const Eigen::VectorXd& x) {
	Eigen::VectorXd grad;
  Eigen::SparseMatrix<double> hessian(x.rows(),x.rows());
  std::vector<Eigen::Triplet<double> > IJV;

  IJV.reserve(quadTop.E.rows()*36);
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0;

  #pragma clang loop vectorize(enable)
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

		double t2 = pxf_x*4.0;
		double t3 = pxf_y*4.0;
		double t4 = pxf_z*4.0;
		double t5 = p0_x*4.0;
		double t10 = p0_x*2.0;
		double t11 = pxf_x*2.0;
		double t6 = t10-t11;
		double t7 = p0_x-pxf_x;
		double t8 = p0_y-pxf_y;
		double t9 = p0_z-pxf_z;
		double t12 = t6*t6;
		double t13 = t12*2.0;
		double t14 = t7*t7;
		double t15 = t14*4.0;
		double t16 = t8*t8;
		double t17 = t16*4.0;
		double t18 = t9*t9;
		double t19 = t18*4.0;
		double t20 = p0_y*2.0;
		double t21 = pxf_y*2.0;
		double t22 = t20-t21;
		double t23 = t6*t22*2.0;
		double t24 = p0_z*2.0;
		double t25 = pxf_z*2.0;
		double t26 = t24-t25;
		double t27 = t6*t26*2.0;
		double t28 = p0_y*4.0;
		double t29 = l0*4.0;
		double t30 = t22*t22;
		double t31 = t30*2.0;
		double t32 = t22*t26*2.0;
		double t33 = p0_z*4.0;
		double t34 = t26*t26;
		double t35 = t34*2.0;
		double t36 = -t2+t5;
		double t37 = -t13-t15-t17-t19+t29;
		double t38 = -t3+t28;
		double t39 = -t15-t17-t19+t29-t31;
		double t40 = t15+t17+t19-t29+t31;
		double t41 = -t4+t33;
		double t42 = -t15-t17-t19+t29-t35;
		double t43 = t15+t17+t19-t29+t35;

		// order is p0_x, p0_y, p0_z, pxf_x, pxf_y, pxf_z

		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_0_i, l0*-4.0+t13+t15+t17+t19));
		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_0_i+vnum, t23));
		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_0_i+2*vnum, t27));
		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_xf_i, t37));
		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_xf_i+vnum, -t23));
		IJV.push_back(Eigen::Triplet<double>(p_0_i,p_xf_i+2*vnum, -t27));

		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_0_i, t23));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t40));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_0_i+2*vnum, t32));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_xf_i, -t23));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t39));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+2*vnum, -t32));

		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i, t27));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+vnum, t32));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t43));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i, -t27));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+vnum, -t32));
		IJV.push_back(Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t42));

		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_0_i, t37));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_0_i+vnum, -t23));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_0_i+2*vnum, -t27));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_xf_i, t13+t15+t17+t19-t29));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_xf_i+vnum, t23));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i,p_xf_i+2*vnum, t27));

		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_0_i, -t23));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t39));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+2*vnum, -t32));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i, t23));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t40));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+2*vnum, t32));

		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i, -t27));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+vnum, -t32));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t42));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i, t27));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+vnum, t32));
		IJV.push_back(Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t43));

        h_cnt++;
  }
  hessian.setFromTriplets(IJV.begin(),IJV.end());
  return hessian;
}