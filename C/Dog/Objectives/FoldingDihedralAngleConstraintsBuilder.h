#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

// For now the class assumes the initialization is flat (so it linearly interpolates dihedral angles from 0 to that)
class FoldingDihedralAngleConstraintsBuilder {
public:
	FoldingDihedralAngleConstraintsBuilder(const Dog& dog, std::vector<double>& destination_dihedral_angles, const double& timestep);
	
	void add_constraint(const EdgePoint& ep, double dihedral_angle);
	void get_edge_angle_pairs(std::vector<std::pair<Edge,Edge>>& out_edge_angle_pairs) {out_edge_angle_pairs = edge_angle_pairs;}
	// convert edge angles to dihedral
	void get_edge_angle_constraints(std::vector<double>& edge_cos_angles);
private:
	void find_prev_next_edge_points(const EdgePoint& ep, EdgePoint& prev_ep, EdgePoint& next_ep);

	const Dog& dog;
	std::vector<double>& destination_dihedral_angles; 
	// between 0 and 1
	const double& timestep;
	std::vector<std::pair<Edge,Edge>> edge_angle_pairs;
	std::vector<double> tangent_angles; // angles between the tangent of the curve and the DOG grid
};
