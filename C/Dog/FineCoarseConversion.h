#pragma once

#include <Eigen/Dense>
#include <vector>

#include "igl/serialize.h"
#include "Dog.h"

using namespace std;

// these would probably be nicer in an enum, but this is how it started...
const int UNCHECKED = -1;
const int FINEONLY_I = -2;//is between two coarse points
const int UNDECIDED = -3;
const int WILLBECHECKED = -4;//used during construction
const int FINEONLY_X = -5;//has only fineonly neighbours
const int COARSEONLY = -6;//only appears in coarse mesh
const int LINK = -7;//is in fine and coarse
const int CURVEONLY_LINK = -8;//is in curve and coarse
const int CURVEONLY_I = -9;//only in curve submesh, between two link vertices

class FineCoarseConversion  : public igl::Serializable {
public:

	FineCoarseConversion(){};
	FineCoarseConversion(const Dog& fine_dog, const Dog& coarse_dog);

	//Get vertex indices in other mesh (negative value if it isn't there)
	int fine_to_coarse(int fine) const {return ftc[fine];}
	int coarse_to_fine(int coarse) const {return ctf[coarse];}

	//Maps certain fine-only points to an edge in the coarse Dog
	int fine_to_coarse_edge(int fine) const {return ftc_edge[fine];}

	int coarse_to_fine_curve(int curve_idx, int ep_idx) const {return ctf_curve[curve_idx][ep_idx];}
	//const vector< vector< vector<double> > >& getCurveOffsets() const {return ctf_curve_offsets;}
	Eigen::MatrixXd getCoarseCurveCoords(const Dog& coarse_dog, int curve_idx) const;
	Eigen::MatrixXd getInterpolatedCurveCoords(const Dog& fine_dog, const Dog& coarse_dog, int curve_idx) const;

	// Computes coarse_V given fine_V, extrapolating the coarseonly vertices
	Eigen::MatrixXd coarsen(const Eigen::MatrixXd& fine_V) const;

	void print() const;

	void InitSerialization() {
		Add(ftc, std::string("fine_to_coarse"));
		Add(ctf, std::string("coarse_to_fine"));
    Add(ftc_edge, std::string("fine_to_coarse_edge"));
    Add(ftc_curve, std::string("fine_to_coarse_stitched_curves"));
    Add(ctf_curve, std::string("coarse_to_fine_stitched_curves"));
    Add(ctf_curve_offsets, std::string("coarse_to_fine_curve_offsets"));
		Add(ftc_update_vertices, std::string("ftc_update_vertices"));
		Add(ftc_update_weights, std::string("ftc_update_weights"));
  }

private:

	//ftc = Fine-To-Coarse
	vector<int> ftc;
	vector<int> ctf;
	vector<int> ftc_edge;

	vector< vector<int> > ftc_curve;
	vector< vector<int> > ctf_curve;
	vector< vector< vector<double> > > ctf_curve_offsets;

	//Data for fine-to-coarse update at coarse-only vertices
	vector< vector<int> > ftc_update_vertices;
	vector< vector<double> > ftc_update_weights;

};
