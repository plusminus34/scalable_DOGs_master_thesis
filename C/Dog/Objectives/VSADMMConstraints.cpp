#include "VSADMMConstraints.h"

using namespace std;

VSADMMConstraints::VSADMMConstraints(): A(nullptr), z(nullptr),
		lambda(nullptr), rho(0) {
	const_n = 0;
	IJV.resize(0);
	ready = false;
}

VSADMMConstraints::VSADMMConstraints(const Eigen::SparseMatrix<double> &A_i,
	  const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i):
 		A(&A_i), z(&z_i), lambda(&lambda_i), rho(rho_i){
	if(A == nullptr) {ready = false;}
	else {
		const_n = A->rows();
		IJV.resize(A->nonZeros());
		ready = true;
	}
}

void VSADMMConstraints::initialize() {
	if (ready) return;
	const_n = A->rows();
	IJV.resize(A->nonZeros());
	ready = true;
}
void VSADMMConstraints::set_pointers(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i){
	A = &A_i;
	z = &z_i;
	lambda = &lambda_i;
	rho = rho_i;
	initialize();
}

Eigen::VectorXd VSADMMConstraints::Vals(const Eigen::VectorXd& x) const {
	return (*A) * x - *z - *lambda / rho;
}

void VSADMMConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {
	if (!ready) return;
	int count = 0;
	for (int i=0; i < A->outerSize(); ++i) {
	  for (Eigen::SparseMatrix<double>::InnerIterator it(*A,i); it; ++it) {
			IJV[count++] = Eigen::Triplet<double>(it.row(), it.col(), it.value());
	  }
	}
}
