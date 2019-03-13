#include "MVFoldingDihedralAngleConstraintsBuilder.h"

using namespace std;

MVFoldingDihedralAngleConstraintsBuilder::MVFoldingDihedralAngleConstraintsBuilder(const Dog& dog, const double& timestep) : 
				dog(dog), timestep(timestep) {
	// Todo: find the initial tangent angles so that we could convert angles to angles
}

void MVFoldingDihedralAngleConstraintsBuilder::add_constraint(const EdgePoint& ep, double dihedral_angle, bool is_mountain_fold) {
	// Get the two tangent edges
	int v1,v2,w1,w2;
	dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);
	Edge e1(v1,v2),e2(w1,w2);

	// Get the next and prev curved fold points
	EdgePoint ep_b,ep_f; find_prev_next_edge_points(ep,ep_b,ep_f);
	Eigen::RowVector3d curve_t1 = ep.getPositionInMesh(dog.getV())-ep_b.getPositionInMesh(dog.getV());
	Eigen::RowVector3d curve_t2 = ep_f.getPositionInMesh(dog.getV())-ep.getPositionInMesh(dog.getV());
	// Take the tangent of the osculating circle passing through these 3 points
	Eigen::RowVector3d T = (curve_t1*curve_t2.norm()+curve_t2*curve_t1.norm()).normalized();

	Eigen::RowVector3d grid_tangent_direction = (dog.getV().row(v1)-dog.getV().row(v2)).normalized();
	//cout << "grid_tangent_direction = " << grid_tangent_direction << endl;
	//cout << "T = " << T << endl;
	// Doesn't matter if we take alpha or pi-alpha as its symmetric and give the same values
	double curve_tangents_angle = acos(T.dot((grid_tangent_direction).normalized()));

	// TODO switch v1,v2 and w1,w2 according to mountain/valley fold
	double t = ep.t;

	// Opposite order for v1,v2 and w1,w2 to make sure there's a fold
	 // note the 't' here is important for the osculating plane, so if changing M/V make sure we change t as well!
	// but 't' is defined on v1,v2 pair (so can change the second set of vertices without it)
	MVTangentCreaseFold mvFold1(v1,v2,w2,w1,t,ep_b,ep_f),mvFold2(w2,w1,v1,v2,ep.t,ep_b,ep_f);
	mvTangentCreaseFolds.push_back(mvFold1);// mvTangentCreaseFolds.push_back(mvFold2);
	is_mountain.push_back(is_mountain_fold); //is_mountain.push_back(is_mountain_fold); // push twice, one for each side of the fold
	tangent_angles.push_back(curve_tangents_angle); //tangent_angles.push_back(curve_tangents_angle); // twice
	destination_dihedral_angles.push_back(dihedral_angle);// destination_dihedral_angles.push_back(dihedral_angle); // twice
	//cout << "Added edge point with v1 = " << v1 << " v2 = " << v2 << " w1 = " << w1 << " w2 = " << w2 << endl;
	//cout << "curve_tangents_angle = " << curve_tangents_angle << endl;

}
void MVFoldingDihedralAngleConstraintsBuilder::find_prev_next_edge_points(const EdgePoint& ep, EdgePoint& prev_ep, EdgePoint& next_ep) {
	const DogEdgeStitching& eS = dog.getEdgeStitching();
	for (int c_i = 0; c_i < eS.stitched_curves.size(); c_i++) {
		for (int e_i = 1; e_i < eS.stitched_curves[c_i].size()-1; e_i++) {
			if (eS.stitched_curves[c_i][e_i].edge == ep.edge) {
				prev_ep = eS.stitched_curves[c_i][e_i-1];
				next_ep = eS.stitched_curves[c_i][e_i+1];
				return;
			}
		}
	}
	std::cout << "Error: could not find prev and next edge points!" << std::endl;
	exit(1); // elegant
}

void MVFoldingDihedralAngleConstraintsBuilder::get_edge_angle_constraints(std::vector<double>& edge_cos_angles) {
	edge_cos_angles.resize(destination_dihedral_angles.size());
	for (int i = 0; i < destination_dihedral_angles.size(); i++) {
		double alpha = tangent_angles[i];
		double dest_dihedral_angle = destination_dihedral_angles[i];
		double const_dihedral = timestep*dest_dihedral_angle;
		std::cout << "const_dihedral = " << const_dihedral << endl;
		// This is half because it is measured from tangent to the bisector
		edge_cos_angles[i] = const_dihedral;//(pow(cos(alpha),2)+pow(sin(alpha),2)*cos(const_dihedral));
		std::cout << "edge_cos_angles[i] = " << edge_cos_angles[i] << std::endl;
		//exit(1);

		//cout << "timestep = " << timestep << endl;
		//cout << "Getting dihedral angle of  " << const_dihedral << " by getting a tangent angle of " << acos(edge_cos_angles[i]) << endl;
	}
}