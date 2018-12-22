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
}

double LaplacianSimilarity::obj(const Eigen::VectorXd& x) const {
	return (L*x-lRef).squaredNorm();
}

Eigen::VectorXd LaplacianSimilarity::grad(const Eigen::VectorXd& x) const {
	Eigen::VectorXd lap_diff = L*x-lRef;
	return 2*(L.transpose()*lap_diff);
}
/*
virtual std::vector<Eigen::Triplet<double> > LaplacianSimilarity::hessianIJV(const Eigen::VectorXd& x) const {
	return L_hessian_IJV;
}
*/
/*
std::vector<Eigen::Triplet<double>> LaplacianSimilarity::to_triplets(Eigen::SparseMatrix<double> & M) {
	std::vector<Eigen::Triplet<double>> v;
    for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    return v;
}
*/


const Eigen::SparseMatrix<double>& LaplacianSimilarity::hessian(const Eigen::VectorXd& x) {
	cachedH = 2*L.transpose();
	return cachedH;
}