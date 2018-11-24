#include "FoldingAnglePositionalConstraintsBuilder.h"

#include "../../GeometryPrimitives/LinearTransformations.h"

using namespace std;
using namespace Eigen;

FoldingAnglePositionalConstraintsBuilder::FoldingAnglePositionalConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS) {
	alpha = 0;
	const_n = 9; // Constrain 4 vertices = 12 vars

	//int c_i = 0.5*eS.edge_const_1.size()/4-3; // TODO: This logic should be inside the constraints builder..
	int c_i = eS.edge_const_1.size()/2;

	edge1 = eS.edge_const_1[c_i]; edge2 = eS.edge_const_2[c_i]; edge_t_coordinate = eS.edge_coordinates[c_i];
	edge1_p1 = V.row(edge1.v1);
	edge1_p2 = V.row(edge1.v2);
	edge2_p1 = V.row(edge2.v1);
	edge2_p2 = V.row(edge2.v2);

	// Rotation center (the point)

	center = EdgePoint(edge1,edge_t_coordinate).getPositionInMesh(V);//edge_t_coordinate * edge1_p1 + (1-edge_t_coordinate) * edge1_p2;
	
	//axis = getRotAxis(V, es, center);
	setRotAxis(V,eS,center);
	std::cout << "axis = " << axis << std::endl;
	
	/*
	Eigen::RowVector3d xAxis(1.,0,0), yAxis(0,1.,0);
	double eps = 1e-12;
	if (abs((edge1_p1-edge1_p2).dot(xAxis)) < eps) {
		axis = xAxis;
	} else {
		axis = yAxis;
	}
	std::cout << "axis = " << axis << std::endl;
	*/
	
	int vnum = V.rows();
	b.resize(const_n);
	//b << edge1.v1,edge1.v1+vnum,edge1.v1+2*vnum, edge1.v2,edge1.v2+vnum,edge1.v2+2*vnum,edge2.v1,edge2.v1+vnum,edge2.v1+2*vnum;//,edge2.v2,edge2.v2+vnum,edge2.v2+2*vnum;
	//b << edge1.v1,edge1.v1+vnum,edge1.v1+2*vnum,edge2.v1,edge2.v1+vnum,edge2.v1+2*vnum,edge2.v2,edge2.v2+vnum,edge2.v2+2*vnum;
	b << edge1.v1,edge2.v1,edge2.v2,
		edge1.v1+vnum,edge2.v1+vnum,edge2.v2+vnum,
		edge1.v1+2*vnum,edge2.v1+2*vnum,edge2.v2+2*vnum;
}

void FoldingAnglePositionalConstraintsBuilder::setRotAxis(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
				const Eigen::RowVector3d& rotCenter) {

	// Find a point a long a curve whose vertex is rotCenter
	for (auto curve: eS.stitched_curves) {
		for (int i = 1; i < curve.size()-1; i++) {
			if (curve[i].getPositionInMesh(V) == rotCenter) {
				std::cout << "Found at i = " << i << std::endl;
				auto pxb = curve[i-1].getPositionInMesh(V), p0 = rotCenter, pxf = curve[i+1].getPositionInMesh(V);
				auto e1 = pxf-p0, e2 = p0-pxb;

				// If the curve is locally straight rotate the edge in the plane
				if (e1.cross(e2).squaredNorm() < 1e-12) {
					// Rotate 90 deg in the plane
					//axis = RowVector3d(e1(1),-e1(0),0.).normalized();
					axis = e1.normalized();
				} else {
					auto tangent = (e1.squaredNorm()*e2 + e2.squaredNorm()*e1).normalized();
					axis = tangent.normalized();
					//auto circle_normal = (e1-e1.dot(tangent)*tangent).normalized();
					//axis = circle_normal;	
				}
				return;
			}
		}
	}
}

void FoldingAnglePositionalConstraintsBuilder::get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) {
	// Rotate edge 1 pt1 around the center (constraint point) in 'axis' with angle alpha
	Eigen::RowVector3d e1_p1_new_loc = rotate_vec(edge1_p1, center, axis, alpha);

	bc_out.resize(const_n);
	bc_out << e1_p1_new_loc(0),edge2_p1(0),edge2_p2(0),
			e1_p1_new_loc(1),edge2_p1(1),edge2_p2(1),
			e1_p1_new_loc(2),edge2_p1(2),edge2_p2(2);

		//bc_out << e1_p1_new_loc(0),e1_p1_new_loc(1),e1_p1_new_loc(2),edge2_p1(0),edge2_p1(1),edge2_p1(2),edge2_p2(0),edge2_p2(1),edge2_p2(2);
	b_out = b;
}
/*
Eigen::VectorXd FoldingAngleConstraints::Vals(const Eigen::VectorXd& x) const {
	Eigen::SparseMatrix<double> Jacobian(const_n, x.rows());
	

	Eigen::VectorXd bc(const_n);
	//bc << e1_p1_new_loc(0),e1_p1_new_loc(1),e1_p1_new_loc(2),edge1_p2(0),edge1_p2(1),edge1_p2(2),edge2_p1(0),edge2_p1(1),edge2_p1(2);//,edge2_p2(0),edge2_p2(1),edge2_p2(2);
	bc << e1_p1_new_loc(0),e1_p1_new_loc(1),e1_p1_new_loc(2),edge2_p1(0),edge2_p1(1),edge2_p1(2),edge2_p2(0),edge2_p2(1),edge2_p2(2);

	Eigen::VectorXd current_vals(const_n);
	for (int i = 0; i < b.rows();i++) {current_vals(i) = x(b(i));}
	auto const_vals = current_vals-bc;
	//cout << "const_vals = " << const_vals << endl;
	return const_vals;
}

std::vector<Eigen::Triplet<double> > FoldingAngleConstraints::JacobianIJV(const Eigen::VectorXd& x) const {
	Eigen::SparseMatrix<double> Jacobian(const_n, x.rows());
	int vnum = x.rows()/3;

	std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz); // All consts are on single coordinates..
	for (int i = 0; i < const_n; i++) {
		IJV.push_back(Eigen::Triplet<double>(i,b(i),1));
	}
	return IJV;
}*/
