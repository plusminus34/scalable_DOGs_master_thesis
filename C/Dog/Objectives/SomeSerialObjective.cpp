#include "SomeSerialObjective.h"

using namespace std;

SomeSerialObjective::SomeSerialObjective(): A(nullptr), lambda(nullptr) {
	IJV.resize(0);
	ready = false;
}

SomeSerialObjective::SomeSerialObjective(const Eigen::SparseMatrix<double> &A_i,
	  const Eigen::VectorXd &lambda_i):
 		A(&A_i), lambda(&lambda_i) {
	if(lambda == nullptr) {ready = false;}
	else {
		IJV.resize(0);
		ready = true;
	}
}

void SomeSerialObjective::initialize() {
	if (ready) return;
	IJV.resize(0);
	ready = true;
}
void SomeSerialObjective::set_pointers(const Eigen::SparseMatrix<double> &A_i,
	 const Eigen::VectorXd &lambda_i) {
	A = &A_i;
	lambda = &lambda_i;
	initialize();
}

double SomeSerialObjective::obj(const Eigen::VectorXd& x) const {
	return - lambda->transpose() * ( (*A) * x );
}

Eigen::VectorXd SomeSerialObjective::grad(const Eigen::VectorXd& x) const {
	return - lambda->transpose() * (*A);
}

void SomeSerialObjective::updateHessianIJV(const Eigen::VectorXd& x) {
}
