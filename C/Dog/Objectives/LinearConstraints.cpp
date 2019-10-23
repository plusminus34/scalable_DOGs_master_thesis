#include "LinearConstraints.h"

using namespace std;

LinearConstraints::LinearConstraints(): A(nullptr), z(nullptr) {
	const_n = 0;
	IJV.resize(0);
	ready = false;
}

LinearConstraints::LinearConstraints(const Eigen::SparseMatrix<double> &A_i,
	  const Eigen::VectorXd &z_i): A(&A_i), z(&z_i) {
	if(A == nullptr) {ready = false;}
	else {
		const_n = A->rows();
		IJV.resize(A->nonZeros());
		ready = true;
	}
}

void LinearConstraints::initialize() {
	if (ready) return;
	const_n = A->rows();
	IJV.resize(A->nonZeros());
	ready = true;
}
void LinearConstraints::set_pointers(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i){
	A = &A_i;
	z = &z_i;
	initialize();
}

Eigen::VectorXd LinearConstraints::Vals(const Eigen::VectorXd& x) const {
	return (*A) * x - *z;//- *lambda / rho;
}

void LinearConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {
	if (!ready) return;
	int count = 0;
	for (int i=0; i < A->outerSize(); ++i) {
	  for (Eigen::SparseMatrix<double>::InnerIterator it(*A,i); it; ++it) {
			IJV[count++] = Eigen::Triplet<double>(it.row(), it.col(), it.value());
	  }
	}
}
