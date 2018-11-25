#include "ModelState.h"

#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/readOBJ.h>

#include "CreasePatterns/SVGReader.h"
#include "CreasePatterns/DogFromCreasePattern.h"

using namespace std;

void ModelState::init_from_mesh(const std::string& mesh_path) {
	std::cout << "Reading mesh " << mesh_path << endl;
	Eigen::MatrixXd V,V_ren; Eigen::MatrixXi F,F_ren;
	igl::readOBJ(mesh_path, V, F_ren);
	
	// We either read a triangle mesh or a quad mesh
	if (F_ren.cols() == 3 ) {
		F = F_to_Fsqr(F_ren);
	} else {
		F = F_ren;F_ren = Fsqr_to_F(F);
	}


	quad_topology(V,F,quadTop);

	// Scale the mesh
	const double edge_l = (V.row(quadTop.bnd_loop[1]) - V.row(quadTop.bnd_loop[0])).norm();
	V *= 1. / edge_l;

	V_ren = V; // no folds
	DogEdgeStitching dogEdgeStitching; // no folds  
	std::vector<pair<int,int>> submesh_min_max; submesh_min_max.push_back(pair<int,int>(0,V.rows()-1)); // only 1 component

	dog = Dog(V,F);
}

void ModelState::init_from_svg(const std::string& svg_path, int x_res, int y_res) {
	CGAL::Bbox_2 bbox;
  	std::vector<Polyline_2> polylines;

	std::cout << "Reading svg " << svg_path << endl;
	read_svg_crease_pattern(svg_path, bbox, polylines);

	CreasePattern creasePattern(bbox, polylines, x_res, y_res);
	dog = dog_from_crease_pattern(creasePattern);
	
	quad_topology(dog.getV(),dog.getF(),quadTop);

  	creasePattern.get_visualization_mesh_and_edges(creasesVisualization.V_arr, creasesVisualization.F_arr, 
											creasesVisualization.faceColors,creasesVisualization.edge_pts1, creasesVisualization.edge_pts2);
  	Eigen::MatrixXd creaseVMesh = dog.getVrendering();
  	double spacing = 3*1.05*CGAL::to_double(bbox.xmax()-bbox.xmin());
  	creaseVMesh.rowwise() += Eigen::RowVector3d(spacing,0,0);
  	Eigen::MatrixXi meshE_i; igl::edges(dog.getFrendering(),meshE_i);
  	Eigen::MatrixXd meshE1; igl::slice(creaseVMesh,meshE_i.col(0),1, creasesVisualization.meshE1);
  	Eigen::MatrixXd meshE2; igl::slice(creaseVMesh,meshE_i.col(1),1, creasesVisualization.meshE2);

	// scale the mesh
	const double edge_l = (dog.getV().row(quadTop.bnd_loop[1]) - dog.getV().row(quadTop.bnd_loop[0])).norm();
	auto scaledV = dog.getV()*1./edge_l;
	dog.update_V(scaledV);
}