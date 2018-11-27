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
	double solve(Dog& dog,
		const PositionalConstraints& posConst,
        const StitchingConstraints& stitchingConstraints,
        /*const EdgePointConstraints& edgePointConstraints, TODO */
        Eigen::MatrixXd& Vout,
        int fixed_mesh_i = -1);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	void procrustes_on_submesh(Dog& dog, int submesh_i, 
								const PositionalConstraints& posConst,
								const EdgePointConstraints& edgePointConstraints);

	void procrustes_on_submesh(Dog& dog, int submesh_i, const PositionalConstraints& posConst);

	PositionalConstraints submesh_positional_constraints_from_mesh_positional_constraints(Dog& dog, int submesh_i, const PositionalConstraints& posConst);
	static int get_best_aligned_submesh(const Dog& dog, const PositionalConstraints& posConst);
	
	int max_iter;
};