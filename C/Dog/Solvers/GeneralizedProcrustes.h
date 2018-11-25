#pragma once

#include "../PositionalConstraints.h"
#include "../EdgePointConstraints.h"

/*
Given n connected components with positional constraints on them, find a rigid motion for each 
disconnected components that match the constraints best.
Also supports edge point constraints.
*/
class GeneralizedProcrustes {
  
public:
	GeneralizedProcrustes(int max_iter): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::MatrixXd& V, const PositionalConstraints& posConst,
        const PositionalConstraints& EdgePointConstraints,
        const PositionalConstraints& posConst, 
        Eigen::MatrixXd& Vout);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	int max_iter;
};