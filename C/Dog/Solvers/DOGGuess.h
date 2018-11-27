#pragma once

#include <igl/arap.h>

#include "../Dog.h"
#include "../Objectives/StitchingConstraints.h"
#include "../../Optimization/PositionalConstraints.h"

class DOGGuess {
public:
	DOGGuess(const Dog& dog, const bool& align_procrustes, const bool& deform_arap);
	void guess(Dog& dog, const PositionalConstraints& postConst, const StitchingConstraints& stitchConst);

	void update_ref(const Eigen::MatrixXd& Vref_i) {Vref = Vref_i;}
private:
	void guessARAP(Dog& dog, const PositionalConstraints& postConst);

	const bool& align_procrustes;
	const bool& deform_arap;

	const Dog& dog_init;
	Eigen::MatrixXd Vref;

	// ARAP related
	Eigen::MatrixXi Ftri; // Triangular mesh for ARAP
	igl::ARAPData arapData;
};