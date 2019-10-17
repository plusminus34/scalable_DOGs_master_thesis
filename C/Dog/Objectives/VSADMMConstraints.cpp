#include "VSADMMConstraints.h"

using namespace std;

VSADMMConstraints::VSADMMConstraints():
 		A(nullptr), z(nullptr), lambda(nullptr), rho(0){			std::cout<<"HIInnn\n";

	const_n = 0;
	IJV.resize(0);
	initd=false;			std::cout<<"HIIaw\n";

}

VSADMMConstraints::VSADMMConstraints(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i):
 		A(&A_i), z(&z_i), lambda(&lambda_i), rho(rho_i){
			std::cout<<"HII\n";
	if(A == nullptr){initd=false;}
	else{
	const_n = A->rows();//z->size();//=x.size() = A.rows()
	IJV.resize(A->nonZeros());
	initd=true;
}	std::cout<<"HIIo\n";

}

void VSADMMConstraints::initnow(){
	if(initd) return;
	const_n = A->rows();//z->size();//=x.size() = A.rows()
	IJV.resize(A->nonZeros());
	initd=true;
}
void VSADMMConstraints::setstuff(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i){
	A = &A_i;
	z = &z_i;
	lambda = &lambda_i;
	rho = rho_i;
	initnow();
}

Eigen::VectorXd VSADMMConstraints::Vals(const Eigen::VectorXd& x) const {
	return (*A) * x - *z - *lambda / rho;//also c/N, but c is 0
}


void VSADMMConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {
	if(!initd) initnow();
	int count = 0;
	for (int i=0; i < A->outerSize(); ++i) {
	  for (Eigen::SparseMatrix<double>::InnerIterator it(*A,i); it; ++it) {
			IJV[count++] = Eigen::Triplet<double>(it.row(), it.col(), it.value());
	  }
	}
}
