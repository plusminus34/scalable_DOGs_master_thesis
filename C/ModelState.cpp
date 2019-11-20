#include "ModelState.h"

#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/readOBJ.h>

#include "CreasePatterns/SVGReader.h"
#include "CreasePatterns/DogFromCreasePattern.h"

using namespace std;

void ModelState::init_from_mesh(const std::string& mesh_path) {
	std::cout << "Reading mesh " << mesh_path << endl;
	Eigen::MatrixXd V; Eigen::MatrixXi F,F_ren;
	igl::readOBJ(mesh_path, V, F_ren);

	// We either read a triangle mesh or a quad mesh
	if (F_ren.cols() == 3 ) {
		F = F_to_Fsqr(F_ren);
	} else {
		F = F_ren;F_ren = Fsqr_to_F(F);
	}
	setup_non_creased_dog(V,F);
}

void ModelState::setup_non_creased_dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
	QuadTopology quadTop; quad_topology(V,F,quadTop);

	// Scale the mesh so that the entire x curve will be of length 20
	const double edge_l = (V.row(quadTop.bnd_loop[1]) - V.row(quadTop.bnd_loop[0])).norm();
	double cur_v_len = edge_l*sqrt(V.rows());
	auto scaled_V =V* (20. / cur_v_len);

	DogEdgeStitching dogEdgeStitching; // no folds
	std::vector<pair<int,int>> submesh_min_max; submesh_min_max.push_back(pair<int,int>(0,V.rows()-1)); // only 1 component

	dog = Dog(scaled_V,F);
}

void ModelState::init_from_planar(int square_h, int square_w) {
	std::cout << "init from planar with  square_h = " << square_h << " square_w = " << square_w << std::endl;
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	get_planar_square_mesh(V, F, square_h, square_w); F = F_to_Fsqr(F);
	setup_non_creased_dog(V,F);
}

void ModelState::init_from_svg(const std::string& svg_path, int x_res, int y_res) {
	CGAL::Bbox_2 bbox;
  	std::vector<Polyline_2> polylines, border_polylines;

	std::cout << "Reading svg " << svg_path << endl;
	read_svg_crease_pattern(svg_path, bbox, polylines, border_polylines, cp_polylines, cp_border_polylines);

	std::cout << "polylines: "<<polylines.size()<<"\n";
	if(polylines.size()>0){
//		std::cout << "polylines[0] is of size "<<polylines[0].rows()<<" x "<<polylines[0].cols()<<"\n";
	}
	std::cout << "border_polylines: "<<border_polylines.size()<<"\n";
	if(border_polylines.size()>0){
//		std::cout << "border_polylines[0] is of size "<<border_polylines[0].rows()<<" x "<<border_polylines[0].cols()<<"\n";
	}

	CreasePattern creasePattern(bbox, polylines, border_polylines, x_res, y_res);
	dog = dog_from_crease_pattern(creasePattern);

  if(x_res%2 == 0 || y_res%2 == 0) cout << "Warning: Resolution should really be odd\n";
	int x_res_coarse = (x_res + 1) / 2;
	int y_res_coarse = (y_res + 1) / 2;
	CreasePattern coarse_creasePattern(bbox, polylines, border_polylines, x_res_coarse, y_res_coarse);
	coarse_dog = dog_from_crease_pattern(coarse_creasePattern);
	conversion = FineCoarseConversion(dog, coarse_dog);

	cp_bb_x_min = bbox.xmin();
	cp_bb_x_max = bbox.xmax();
	cp_bb_y_min = bbox.ymin();
	cp_bb_y_max = bbox.ymax();
	cp_x_res = x_res;
	cp_y_res = y_res;

	/*
  	creasePattern.get_visualization_mesh_and_edges(creasesVisualization.V_arr, creasesVisualization.F_arr,
											creasesVisualization.faceColors,creasesVisualization.edge_pts1, creasesVisualization.edge_pts2);
											*/
	creasePattern.get_visualization_edges(creasesVisualization.edge_pts1, creasesVisualization.edge_pts2);
  	Eigen::MatrixXd creaseVMesh = dog.getVrendering();
  	double spacing = 3*1.05*CGAL::to_double(bbox.xmax()-bbox.xmin());
  	creaseVMesh.rowwise() += Eigen::RowVector3d(spacing,0,0);
  	Eigen::MatrixXi meshE_i; igl::edges(dog.getFrendering(),meshE_i);
  	Eigen::MatrixXd meshE1; igl::slice(creaseVMesh,meshE_i.col(0),1, creasesVisualization.meshE1);
  	Eigen::MatrixXd meshE2; igl::slice(creaseVMesh,meshE_i.col(1),1, creasesVisualization.meshE2);

	// scale the mesh
	const QuadTopology& quadTop = dog.getQuadTopology();
	const double edge_l = (dog.getV().row(quadTop.bnd_loop[1]) - dog.getV().row(quadTop.bnd_loop[0])).norm();
	auto scaledV = dog.getV()*1./edge_l;
	dog.update_V(scaledV);
	auto coarse_scaledV = (coarse_dog.getV()*1./edge_l) *0.5;//coarse scale
	coarse_dog.update_V(coarse_scaledV);
}
