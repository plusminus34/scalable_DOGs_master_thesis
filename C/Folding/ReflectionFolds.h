#pragma once

#include "../Dog/Dog.h"
#include "../QuadMesh/Quad.h"


// Creates a positional constraint for one vertex, to be in the reflection of the current other vertex in the current osculating plane (with the correct edge length)
// v1 should always be the constant vertex. These constraints should be constructed in some order, for instance by BFS on the submesh's connectivitiy graph
struct ReflectionFold {
	ReflectionFold(const Dog& dog, int curve_idx, int edge_idx, std::vector<bool>& submeshes_set);

	void get_constraint_indices(const Dog& dog, Eigen::VectorXi& b, EdgePoint& edgePoint) const;
	void get_constraint_coords(const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords) const;

	EdgePoint ep, ep_b, ep_f;
	int v1,v2; double len2;
};

class ReflectionFoldConstraintsBuilder {
public:
	void add_fold(const Dog& dog, int curve_idx, int edge_idx, std::vector<bool>& submeshes_set);
	int get_folds_num() const {return folds.size();}
	void clear_folds() {folds.clear();}

	void get_folds_constraint_indices(const Dog& dog, Eigen::VectorXi& b) const;
	void get_folds_constraint_coords(const Dog& dog, Eigen::VectorXd& bc) const;

	std::vector<ReflectionFold> folds;
};