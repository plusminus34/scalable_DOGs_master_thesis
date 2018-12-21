#pragma once

#include "../Dog.h"
#include "../Objectives/StitchingConstraints.h"
#include "../../Optimization/PositionalConstraints.h"
#include "../../Optimization/EdgePointConstraints.h"

class DOGGuess {
public:
	DOGGuess(const Dog& dog, const bool& align_procrustes);
	void guess(Dog& dog, const PositionalConstraints& posConst, StitchingConstraints& stitchConst,
		EdgePointConstraints& edgePointConstraints);
private:

	const bool& align_procrustes;
};