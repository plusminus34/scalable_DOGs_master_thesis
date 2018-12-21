#include "DOGGuess.h"

#include "GeneralizedProcrustes.h"
#include "../../Optimization/CompositeConstraints.h"

using namespace igl;

DOGGuess::DOGGuess(const Dog& dog, const bool& align_procrustes) : align_procrustes(align_procrustes) {}

void DOGGuess::guess(Dog& dog, const PositionalConstraints& posConst, StitchingConstraints& stitchConst,
		EdgePointConstraints& edgePointConstraints) {
	if (align_procrustes){
		GeneralizedProcrustes genProc; genProc.solve(dog, posConst, stitchConst, edgePointConstraints);	
	}
	dog.update_Vren();
}