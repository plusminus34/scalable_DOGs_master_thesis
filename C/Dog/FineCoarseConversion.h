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

	//Get vertex indices in other mesh (negative value if it isn't there)
	int fine_to_coarse(int fine) const {return ftc(fine);}
	int coarse_to_fine(int coarse) const {return ctf(coarse);}

	//Maps certain fine-only points to an edge in the coarse Dog
	int fine_to_coarse_edge(int fine) const {return ftc_edge[fine];}

	int coarse_to_fine_curve(int curve_idx, int ep_idx);
	vector< vector< vector<double> > >& getCurveOffsets(){return ctf_curve_offsets;}

	void print();

	//This one is deprecated
	Eigen::RowVector3d get_fine_curve_approx(const Eigen::MatrixXd& V, int curve_idx, int ep_idx);

	void InitSerialization() {
		Add(ftc, std::string("fine_to_coarse"));
		Add(ctf, std::string("coarse_to_fine"));
    Add(ftc_curve, std::string("fine_to_coarse_stitched_curves"));
    Add(ctf_curve, std::string("coarse_to_fine_stitched_curves"));
		Add(ctf_curve_vertices, std::string("curve_interpolation_vertices"));
		Add(ctf_curve_weights, std::string("curve_interpolation_weights"));
		Add(ctf_curve_offsets, std::string("curve_offsets"));
  }

private:

	//ftc = Fine-To-Coarse
	Eigen::VectorXi ftc;
	Eigen::VectorXi ctf;
	Eigen::VectorXi ftc_edge;

	vector< vector<int> > ftc_curve;
	vector< vector<int> > ctf_curve;
	vector< Eigen::MatrixXi > ctf_curve_vertices;
	vector< Eigen::MatrixXd > ctf_curve_weights;

	vector< vector< vector<double> > > ctf_curve_offsets;
};
