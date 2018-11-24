#pragma once

#include <igl/arap.h>

#include "../Dog.h"
#include "../../Optimization/PositionalConstraints.h"

class DOGGuess {
public:
	DOGGuess(const Dog& dog);
	void guess(Eigen::MatrixXd& V, const PositionalConstraints& postConst, Eigen::MatrixXd& guess);

	void update_ref(const Eigen::MatrixXd& Vref_i) {Vref = Vref_i;}
private:
	void guessARAP(Eigen::MatrixXd& V, const PositionalConstraints& postConst, Eigen::MatrixXd& guess);
	//void guessProcrustes(Eigen::MatrixXd& V, const PositionalConstraints& postConst);

	const Dog& dog_init;
	Eigen::MatrixXd Vref;

	// ARAP related
	Eigen::MatrixXi Ftri; // Triangular mesh for ARAP
	igl::ARAPData arapData;
};