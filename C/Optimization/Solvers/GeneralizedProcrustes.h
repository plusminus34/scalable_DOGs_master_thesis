#pragma once

#include "../PositionalConstraints.h"
#include "../EdgePointConstraints.h"

/*
Given n connected components with positional constraints on them, find a rigid motion for each 
disconnected components that match the constraints best.
Supports also dog edge stitching constraints, and later maybe edge point constraints.

As it works with propagation, get the "first" mesh as input. 
This one just match a rigid motion for the positional constraints.
Then it "fixes" the values for the linear constraints of nearby values, for instance stitching constraints.
*/
class GeneralizedProcrustes {
  
public:
	GeneralizedProcrustes(int max_iter): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::MatrixXd& V, const PositionalConstraints& posConst,
        const PositionalConstraints& posConst,
        const StitchingConstraints& stitchingConst, 
        Eigen::MatrixXd& Vout);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	int max_iter;
};