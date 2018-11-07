#include "Dog.h"

Dog::Dog(Eigen::MatrixXd V, Eigen::MatrixXi F, DogFoldingConstraints foldingConstraints, Eigen::MatrixXi F_ren) : 
				V(V),F(F),foldingConstraints(foldingConstraints),F_ren(F_ren) {
	// empty
	// TODO: create/update V_ren?
}

Dog::Dog(const Dog& d) : V(d.V),F(d.F),foldingConstraints(d.foldingConstraints),F_ren(d.F_ren){
	// empty
	// TODO: create/update V_ren?
}

void Dog::get_V_ren(const Eigen::MatrixXd& V, const DogFoldingConstraints& fC, Eigen::MatrixXd& V_ren) {
	int consts_num = fC.edge_coordinates.size();
	Eigen::MatrixXd V_folds_polygons(consts_num,3);
	for (int const_i = 0; const_i < consts_num; const_i++) {
		double t = fC.edge_coordinates[const_i];
		V_folds_polygons.row(const_i) = t*V.row(fC.edge_const_1[const_i].first) + (1-t)*V.row(fC.edge_const_1[const_i].second);
	}
	V_ren.resize(V.rows()+consts_num,3);
	V_ren << V,V_folds_polygons;
}