#include "DOGGuess.h"

#include "GeneralizedProcrustes.h"

DOGGuess::DOGGuess(const Dog& dog, const bool& align_procrustes, const bool& deform_arap) : dog_init(dog), align_procrustes(align_procrustes),
															deform_arap(deform_arap) {
	Ftri = Fsqr_to_F(dog.getF());
	arapData.max_iter = 5;	
}

void DOGGuess::guess(Dog& dog, const PositionalConstraints& postConst, const StitchingConstraints& stitchConst) {
	if (align_procrustes){
		GeneralizedProcrustes genProc; genProc.solve(dog, postConst, stitchConst);	
	}
	if (deform_arap) {
		guessARAP(dog, postConst);
	}
}

void DOGGuess::guessARAP(Dog& dog, const PositionalConstraints& postConst) {
	auto b = postConst.getPositionIndices(); auto bc = postConst.getPositionVals();
	Eigen::VectorXi b_V(b.rows()/3); Eigen::MatrixXd bc_V(b_V.rows(),3);
	for (int i =0 ; i < b_V.rows(); i++) {b_V(i) = b(i);}
	vec_to_mat2(bc, bc_V);

	// getPositionIndices,getPositionBals
	igl::arap_precomputation(Vref,Ftri,3,b_V,arapData);
	igl::arap_solve(bc_V,arapData,dog.getVMutable());
}