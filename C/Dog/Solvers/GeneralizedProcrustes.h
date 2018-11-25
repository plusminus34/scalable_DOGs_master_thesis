#pragma once

#include "../../Optimization/PositionalConstraints.h"
#include "../../Optimization/EdgePointConstraints.h"
#include "../Objectives/StitchingConstraints.h"

/*
Given n connected components with positional constraints on them, find a rigid motion for each 
disconnected components that match the constraints best.
Also supports edge point constraints.
*/
class GeneralizedProcrustes {
  
public:
	GeneralizedProcrustes(int max_iter): max_iter(max_iter) {}
	// Fixed mesh is the input of the first submesh that should fix vertices for its neighbouring submeshes
	double solve(Dog& dog, int fixed_mesh_i,
		const PositionalConstraints& posConst,
        const StitchingConstraints& stitchingConstraints,
        /*const PositionalConstraints& EdgePointConstraints, TODO */
        Eigen::MatrixXd& Vout);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	void procrustes_on_submesh(Dog& dog, int fixed_mesh_i, const PositionalConstraints& posConst);
	int max_iter;
};