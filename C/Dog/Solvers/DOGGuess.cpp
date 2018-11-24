#include "DOGGuess.h"

DOGGuess::DOGGuess(const Dog& dog) : dog_init(dog) {
	Ftri = Fsqr_to_F(dog.getF());
	arapData.max_iter = 5;	
}

void DOGGuess::guess(Eigen::MatrixXd& V, const PositionalConstraints& postConst, Eigen::MatrixXd& guess) {
	guessARAP(V, postConst, guess);
}

void DOGGuess::guessARAP(Eigen::MatrixXd& V, const PositionalConstraints& postConst, Eigen::MatrixXd& guess) {
	auto b = postConst.getPositionIndices(); auto bc = postConst.getPositionVals();
	Eigen::VectorXi b_V(b.rows()/3); Eigen::MatrixXd bc_V(b_V.rows(),3);
	for (int i =0 ; i < b_V.rows(); i++) {b_V(i) = b(i);}
	vec_to_mat2(bc, bc_V);

	// getPositionIndices,getPositionBals
	igl::arap_precomputation(Vref,Ftri,3,b_V,arapData);
	igl::arap_solve(bc_V,arapData,guess);
}