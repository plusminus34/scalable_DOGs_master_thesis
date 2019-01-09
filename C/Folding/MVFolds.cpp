#include "MVFolds.h"

using namespace std;

MountainValleyFold::MountainValleyFold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain) : is_mountain(is_mountain) {
	auto eS = dog.getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_idx];	
	ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_b = foldingCurve[edge_idx+1];
	dog.get_2_inner_vertices_from_edge(ep.edge,v1,v2);

	Eigen::RowVector3d edgeCoords = ep.getPositionInMesh(dog.getV());
	len1 = (edgeCoords-dog.getV().row(v1)).norm();
	len2 = (edgeCoords-dog.getV().row(v2)).norm();

}
void MVFoldingConstraintsBuilder::add_fold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain) {
	MountainValleyFold mvF(dog, curve_idx, edge_idx, is_mountain);
	folds.push_back(mvF);
}
/*
struct MountainValleyFold {
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
	void get_folds_constraint_coords(const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords);
private:
	std::vector<MountainValleyFold> folds;
};
*/