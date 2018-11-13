#include "ModelState.h"

#include <igl/readOBJ.h>

using namespace std;

void ModelState::init_from_mesh(const std::string& mesh_path) {
	std::cout << "Reading mesh " << mesh_path << endl;
	Eigen::MatrixXd V,V_ren; Eigen::MatrixXi F,F_ren;
	igl::readOBJ(mesh_path, V, F_ren);
	
	if (F_ren.cols() == 3 ) {
		// We read a triangle mesh
		F = F_to_Fsqr(F_ren);
	} else {
		// We read a quad mesh
		F = F_ren;
		F_ren = Fsqr_to_F(F);
	}
	quad_topology(V,F,quadTop);

	// Scale the mesh
	const double edge_l = (V.row(quadTop.bnd_loop[1]) - V.row(quadTop.bnd_loop[0])).norm();
	V *= 1. / edge_l;

	V_ren = V; // no folds
	DogEdgeStitching dogEdgeStitching; // no folds  

	dog = Dog(V,F,dogEdgeStitching,V_ren,F_ren);
}