#include "SimplifiedBendingObjective.h"

SimplifiedBendingObjective::SimplifiedBendingObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {
	// Number of hessian triplets
	IJV.resize(51*quadTop.stars.rows()/5+27*quadTop.bnd3.rows()/4);
}

double SimplifiedBendingObjective::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  int vnum = x.rows()/3;
    
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double t5 = p0_x*2.0;
		double t2 = pxb_x+pxf_x-t5;
		double t7 = p0_y*2.0;
		double t3 = pxb_y+pxf_y-t7;
		double t9 = p0_z*2.0;
		double t4 = pxb_z+pxf_z-t9;
		double t6 = pyb_x+pyf_x-t5;
		double t8 = pyb_y+pyf_y-t7;
		double t10 = pyb_z+pyf_z-t9;
		e += t2*t2+t3*t3+t4*t4+t6*t6+t8*t8+t10*t10;

  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));
    
	double t2 = p0_x*-2.0+pxb_x+pxf_x;
	double t3 = p0_y*-2.0+pxb_y+pxf_y;
	double t4 = p0_z*-2.0+pxb_z+pxf_z;
  	e += t2*t2+t3*t3+t4*t4;
  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  return e;
}

Eigen::VectorXd SimplifiedBendingObjective::grad(const Eigen::VectorXd& x) const {
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

        double t2 = pxb_x*2.0;
		double t3 = pxf_x*2.0;
		double t11 = p0_x*4.0;
		double t4 = t2+t3-t11;
		double t5 = pxb_y*2.0;
		double t6 = pxf_y*2.0;
		double t12 = p0_y*4.0;
		double t7 = t5+t6-t12;
		double t8 = pxb_z*2.0;
		double t9 = pxf_z*2.0;
		double t13 = p0_z*4.0;
		double t10 = t8+t9-t13;
		double t14 = pyb_x*2.0;
		double t15 = pyf_x*2.0;
		double t16 = -t11+t14+t15;
		double t17 = pyb_y*2.0;
		double t18 = pyf_y*2.0;
		double t19 = -t12+t17+t18;
		double t20 = pyb_z*2.0;
		double t21 = pyf_z*2.0;
		double t22 = -t13+t20+t21;

		grad(p_0_i) += p0_x*1.6E1-pxb_x*4.0-pxf_x*4.0-pyb_x*4.0-pyf_x*4.0;
		grad(p_0_i+v_num) += p0_y*1.6E1-pxb_y*4.0-pxf_y*4.0-pyb_y*4.0-pyf_y*4.0;
		grad(p_0_i+2*v_num) += p0_z*1.6E1-pxb_z*4.0-pxf_z*4.0-pyb_z*4.0-pyf_z*4.0;
		grad(p_xb_i+0) += t4;
		grad(p_xb_i+v_num) += t7;
		grad(p_xb_i+2*v_num) += t10;
		grad(p_xf_i) += t4;
		grad(p_xf_i+v_num) += t7;
		grad(p_xf_i+2*v_num) += t10;
		grad(p_yb_i) += t16;
		grad(p_yb_i+v_num) += t19;
		grad(p_yb_i+2*v_num) += t22;
		grad(p_yf_i) += t16;
		grad(p_yf_i+v_num) += t19;
		grad(p_yf_i+2*v_num) += t22;
    //}
  }


  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	double t2 = pxb_x*2.0;
	double t3 = pxf_x*2.0;
	double t4 = p0_x*-4.0+t2+t3;
	double t5 = pxb_y*2.0;
	double t6 = pxf_y*2.0;
	double t7 = p0_y*-4.0+t5+t6;
	double t8 = pxb_z*2.0;
	double t9 = pxf_z*2.0;
	double t10 = p0_z*-4.0+t8+t9;
	
	grad(p_0_i) += p0_x*8.0-pxb_x*4.0-pxf_x*4.0;
	grad(p_0_i+v_num) += p0_y*8.0-pxb_y*4.0-pxf_y*4.0;
	grad(p_0_i+2*v_num) += p0_z*8.0-pxb_z*4.0-pxf_z*4.0;
	grad(p_xb_i+0) += t4;
	grad(p_xb_i+v_num) += t7;
	grad(p_xb_i+2*v_num) += t10;
	grad(p_xf_i) += t4;
	grad(p_xf_i+v_num) += t7;
	grad(p_xf_i+2*v_num) += t10;

  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  //cout << "new grad.norm() = " << grad.norm() << endl; //exit(1);
  return grad;
}


void SimplifiedBendingObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  // Number of ijv values
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0;
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

		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i, 1.6E1);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xb_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_yb_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_yf_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, 1.6E1);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_yb_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_yf_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, 1.6E1);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yb_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yf_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_yb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i,p_yf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yf_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_yb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i,p_yf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yf_i+2*vnum, 2.0);
   }

   for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
		int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i, 8.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xb_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, 8.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, 8.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, -4.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, 2.0);
		IJV[h_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, 2.0);

	}
}