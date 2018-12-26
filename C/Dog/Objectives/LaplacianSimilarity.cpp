#include "LaplacianSimilarity.h"

#include "../DogLaplacian.h"

LaplacianSimilarity::LaplacianSimilarity(const Dog& dog, const Eigen::VectorXd& x0) {
	Eigen::MatrixXd V_x; vec_to_mat2(x0,V_x); 
	L = DOG_laplacian(V_x,dog.getF());

	/*
	// The hessian is 2*L
	L_hessian_IJV = to_triplets(L); 
	// Now double the values
	for (int i = 0; i < L_hessian_IJV.size(); i++) {
		L_hessian_IJV[i] = Eigen::Triplet<double>(L_hessian_IJV[i].row(),L_hessian_IJV[i].col(),2*L_hessian_IJV[i].value());
	}
	*/
	lRef = L*x0;
	cachedH = 2*L.transpose();
	IJV = to_triplets(cachedH);
}

double LaplacianSimilarity::obj(const Eigen::VectorXd& x) const {
	return (L*x-lRef).squaredNorm();
}

Eigen::VectorXd LaplacianSimilarity::grad(const Eigen::VectorXd& x) const {
	Eigen::VectorXd lap_diff = L*x-lRef;
	return 2*(L.transpose()*lap_diff);
}

// TODO: This completely ignores the hessian of the constraints! (works perfectly for linear constraints such as positions though)
void LaplacianSimilarity::updateHessianIJV(const Eigen::VectorXd& x) {
	// Could write it directly maybe, or have an eddificent A*A' at least..
	// Or preallocate the IJV (second time)
	IJV = to_triplets(cachedH);
}


const Eigen::SparseMatrix<double>& LaplacianSimilarity::hessian(const Eigen::VectorXd& x) {
	return cachedH;
}