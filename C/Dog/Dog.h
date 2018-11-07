#pragma once

#include <Eigen/Dense>
#include <vector>

struct DogFoldingConstraints {
	std::vector<std::pair<int,int>> edge_const_1, edge_const_2;
	std::vector<double> edge_coordinates;
};

class Dog {
public:
	Dog(Eigen::MatrixXd V, Eigen::MatrixXi F, DogFoldingConstraints foldingConstraints, Eigen::MatrixXi F_ren);
	
private:
	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	// The initial rendered mesh
	Eigen::MatrixXi F_ren;

	// Folding constraints
	DogFoldingConstraints foldingConstraints;
};