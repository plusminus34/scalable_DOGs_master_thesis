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
	Dog(const Dog& dog);

	static void get_V_ren(const Eigen::MatrixXd& V, const DogFoldingConstraints& foldingConstraints, Eigen::MatrixXd& V_ren);
	
private:
	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	// The initial rendered (triangular) mesh
	Eigen::MatrixXi F_ren;

	// Folding constraints
	DogFoldingConstraints foldingConstraints;
};