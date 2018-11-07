#include "Dog.h"

Dog::Dog(Eigen::MatrixXd V, Eigen::MatrixXi F, DogFoldingConstraints foldingConstraints, Eigen::MatrixXi F_ren) : 
				V(V),F(F),foldingConstraints(foldingConstraints),F_ren(F_ren) {
	// empty
	// TODO: create/update V_ren?
}
