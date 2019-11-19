#pragma once

#include <Eigen/Dense>
#include <vector>

#include "igl/serialize.h"
#include "Dog.h"

using namespace std;

class FineCoarseConversion  : public igl::Serializable {
public:

	FineCoarseConversion(){};
	FineCoarseConversion(const Dog& fine_dog, const Dog& coarse_dog);

	//Get vertex indices in other mesh (-1 if it isn't there)
	int fine_to_coarse(int fine);
	int coarse_to_fine(int coarse);

	int coarse_to_fine_curve(int curve_idx, int ep_idx);

	Eigen::RowVector3d get_fine_approx(const Eigen::MatrixXd& V, int fine_v);
	Eigen::RowVector3d get_fine_curve_approx(const Eigen::MatrixXd& V, int curve_idx, int ep_idx);

	void InitSerialization() {
		Add(ftc, std::string("fine_to_coarse"));
		Add(ctf, std::string("coarse_to_fine"));
    Add(ftc_curve, std::string("fine_to_coarse_stitched_curves"));
    Add(ctf_curve, std::string("coarse_to_fine_stitched_curves"));
		Add(ctf_curve_vertices, std::string("curve_interpolation_vertices"));
		Add(ctf_curve_weights, std::string("curve_interpolation_weights"));
  }

private:

	//ftc = Fine-To-Coarse
	Eigen::VectorXi ftc;
	Eigen::VectorXi ctf;

	vector< vector<int> > ftc_curve;
	vector< vector<int> > ctf_curve;
	vector< Eigen::MatrixXi > ctf_curve_vertices;
	vector< Eigen::MatrixXd > ctf_curve_weights;
};
