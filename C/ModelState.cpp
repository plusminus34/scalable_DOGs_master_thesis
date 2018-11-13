#include "ModelState.h"

#include <igl/readOBJ.h>

#include "CreasePatterns/SVGReader.h"
#include "CreasePatterns/DogFromCreasePattern.h"

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

void ModelState::init_from_svg(const std::string& svg_path, int x_res, int y_res) {
	CGAL::Bbox_2 bbox;
  	std::vector<Polyline_2> polylines;

	std::cout << "Reading svg " << svg_path << endl;
	read_svg_crease_pattern(svg_path, bbox, polylines);

	CreasePattern creasePattern(bbox, polylines, x_res, y_res);
	dog = dog_from_crease_pattern(creasePattern);
	
	quad_topology(dog.getV(),dog.getF(),quadTop);

	// scale the mesh
	const double edge_l = (dog.getV().row(quadTop.bnd_loop[1]) - dog.getV().row(quadTop.bnd_loop[0])).norm();
	auto scaledV = dog.getV()*1./edge_l;
	dog.update_V(scaledV);

	creasePattern.get_visualization_mesh_and_edges(creasesVisualization.V_arr, creasesVisualization.F_arr, 
											creasesVisualization.faceColors,creasesVisualization.edge_pts1, creasesVisualization.edge_pts2);
}