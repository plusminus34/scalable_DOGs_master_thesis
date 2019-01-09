#pragma once

#include "../Dog/Dog.h"
#include "../QuadMesh/Quad.h"


// is_mountain might need to be checked using the normal z or minus z at the beginning (or deviation from it)
// It will probably be enough to check the orientation of txprincipal_n in the constructor
// At that case if it is flipped just change the is_mountain to be opposite than the input (probably the easiest solution)
struct MountainValleyFold {
	MountainValleyFold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain);

	void get_constraint_indices(const Dog& dog, Eigen::VectorXi& b, EdgePoint& edgePoint);
	void get_constraint_coords(double folding_angle, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords);

	bool is_mountain;
	EdgePoint ep, ep_b, ep_f;
	int v1,v2; double len1, len2;
	double curve_tangents_angle;
};

class MVFoldingConstraintsBuilder {
public:
	void add_fold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain);
	int get_folds_num() {return folds.size();}
	void clear_folds() {folds.clear();}

	void get_folds_constraint_indices(const Dog& dog, Eigen::VectorXi& b, std::vector<EdgePoint>& edgePoints);
	void get_folds_constraint_coords(std::vector<double>& folding_angles, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords);
	// Same angle for all folds
	void get_folds_constraint_coords(double folding_angle, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords);
private:
	std::vector<MountainValleyFold> folds;
};