#include "ProximalObjective.h"

using namespace std;

ProximalObjective::ProximalObjective(): P(nullptr), x_old(nullptr) {
	IJV.resize(0);
	ready = false;
}

ProximalObjective::ProximalObjective(const Eigen::SparseMatrix<double> &P_i,
	  const Eigen::VectorXd &x_old_i):
 		P(&P_i), x_old(&x_old_i) {
	if(P == nullptr) {ready = false;}
	else {
		IJV.resize(P->nonZeros());
		ready = true;
	}
}

void ProximalObjective::initialize() {
	if (ready) return;
	IJV.resize(P->nonZeros());
	ready = true;
}
void ProximalObjective::set_pointers(const Eigen::SparseMatrix<double> &P_i,
	 const Eigen::VectorXd &x_old_i) {
	P = &P_i;
	x_old = &x_old_i;
	initialize();
}

double ProximalObjective::obj(const Eigen::VectorXd& x) const {
	Eigen::VectorXd x_diff = x - *x_old;
	return 0.5 * x_diff.transpose() * (*P) * x_diff;
}

Eigen::VectorXd ProximalObjective::grad(const Eigen::VectorXd& x) const {
	return (*P) * (x - *x_old);
}

void ProximalObjective::updateHessianIJV(const Eigen::VectorXd& x) {
	if (!ready) return;
	int count = 0;
	for (int i=0; i < P->outerSize(); ++i) {
	  for (Eigen::SparseMatrix<double>::InnerIterator it(*P,i); it; ++it) {
			IJV[count++] = Eigen::Triplet<double>(it.row(), it.col(), it.value());
	  }
	}
}
