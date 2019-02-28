#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class FoldingDihedralAngleConstraintsBuilder {
public:
	FoldingDihedralAngleConstraintsBuilder(const Dog& dog, std::vector<double>& destination_dihedral_angles, const double& timestep);
	
	void get_edge_angle_pairs(std::vector<std::pair<Edge,Edge>>& out_edge_angle_pairs) {out_edge_angle_pairs = edge_angle_pairs;}
	// convert edge angles to dihedral
	void get_edge_angle_constraints(std::vector<double>& edge_cos_angles);
private:
	const Dog& dog;
	std::vector<double>& destination_dihedral_angles; 
	// between 0 and 1
	const double& timestep;
	std::vector<std::pair<Edge,Edge>> edge_angle_pairs;
	
};
