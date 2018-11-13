#include "LaplacianSimilarity.h"

#include "../DogLaplacian.h"

LaplacianSimilarity::LaplacianSimilarity(const Dog& dog, const Eigen::VectorXd& x0) {
	Eigen::MatrixXd V_x; vec_to_mat2(x0,V_x); 
	L = DOG_laplacian(V_x,dog.getF());
	lRef = L*x0;
}

double LaplacianSimilarity::obj(const Eigen::VectorXd& x) const {
	return (L*x-lRef).squaredNorm();
}

Eigen::VectorXd LaplacianSimilarity::grad(const Eigen::VectorXd& x) const {
	Eigen::VectorXd lap_diff = L*x-lRef;
	return 2*(L.transpose()*lap_diff);
}
Eigen::SparseMatrix<double> LaplacianSimilarity::hessian(const Eigen::VectorXd& x) const {
	return 2*L.transpose();
}
