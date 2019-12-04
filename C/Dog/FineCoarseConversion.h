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
	int fine_to_coarse_edge(int fine) const {return ftc_edge(fine);}

	int coarse_to_fine_curve(int curve_idx, int ep_idx) const;
	const vector< vector< vector<double> > >& getCurveOffsets() const {return ctf_curve_offsets;}
	Eigen::MatrixXd getInterpolatedCurveCoords(const Dog& fine_dog, const Dog& coarse_dog, int curve_idx) const;

	vector<int> getFineNeighboursOfCoarseOnly(int coarse) const;

	void print() const;

	void InitSerialization() {
		Add(ftc, std::string("fine_to_coarse"));
		Add(ctf, std::string("coarse_to_fine"));
    Add(ftc_edge, std::string("fine_to_coarse_edge"));
    Add(ftc_curve, std::string("fine_to_coarse_stitched_curves"));
    Add(ctf_curve, std::string("coarse_to_fine_stitched_curves"));

		Add(coarseonly_adjacent_links, std::string("CO_adjacent_links"));
		Add(coarseonly_adjacent_fineonly, std::string("CO_adjacent_FOIs"));
  }

private:

	//ftc = Fine-To-Coarse
	Eigen::VectorXi ftc;
	Eigen::VectorXi ctf;
	Eigen::VectorXi ftc_edge;

	vector< vector<int> > ftc_curve;
	vector< vector<int> > ctf_curve;
	vector< vector< vector<double> > > ctf_curve_offsets;

	vector< vector<int> > coarseonly_adjacent_links;//link indices from fine mesh
	vector< vector<int> > coarseonly_adjacent_fineonly;
};
