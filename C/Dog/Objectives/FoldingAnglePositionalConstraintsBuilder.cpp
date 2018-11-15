#include "FoldingAnglePositionalConstraintsBuilder.h"

#include <Eigen/Geometry> 

using namespace std;

FoldingAnglePositionalConstraintsBuilder::FoldingAnglePositionalConstraintsBuilder(const Eigen::MatrixXd& V, Edge i_edge1, Edge i_edge2, 
										double edge_t_coordinate) : edge_t_coordinate(edge_t_coordinate) {
	alpha = 0;
	const_n = 9; // Constrain 4 vertices = 12 vars
	edge1 = i_edge1; edge2 = i_edge2;
	edge1_p1 = V.row(edge1.v1);
	edge1_p2 = V.row(edge1.v2);
	edge2_p1 = V.row(edge2.v1);
	edge2_p2 = V.row(edge2.v2);
	
	Eigen::RowVector3d xAxis(1.,0,0), yAxis(0,1.,0);
	double eps = 1e-12;
	if (abs((edge1_p1-edge1_p2).dot(xAxis)) < eps) {
		axis = xAxis;
	} else {
		axis = yAxis;
	}
	center = edge_t_coordinate * edge1_p1 + (1-edge_t_coordinate) * edge1_p2;
	int vnum = V.rows();
	b.resize(const_n);
	//b << edge1.v1,edge1.v1+vnum,edge1.v1+2*vnum, edge1.v2,edge1.v2+vnum,edge1.v2+2*vnum,edge2.v1,edge2.v1+vnum,edge2.v1+2*vnum;//,edge2.v2,edge2.v2+vnum,edge2.v2+2*vnum;
	b << edge1.v1,edge1.v1+vnum,edge1.v1+2*vnum,edge2.v1,edge2.v1+vnum,edge2.v1+2*vnum,edge2.v2,edge2.v2+vnum,edge2.v2+2*vnum;
}

void FoldingAnglePositionalConstraintsBuilder::get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) {
	// Rotate edge 1 pt1 around the center (constraint point) in 'axis' with angle alpha
	Eigen::RowVector3d e1_p1_new_loc = rotate_vec(edge1_p1, center, axis, alpha);

	bc_out.resize(const_n);
	bc_out << e1_p1_new_loc(0),e1_p1_new_loc(1),e1_p1_new_loc(2),edge2_p1(0),edge2_p1(1),edge2_p1(2),edge2_p2(0),edge2_p2(1),edge2_p2(2);
	b_out = b;

	//Eigen::VectorXd current_vals(const_n);
	//for (int i = 0; i < b.rows();i++) {current_vals(i) = x(b(i));}
	std::cout << "bc = " << bc_out << endl;
	//std::cout << "current vals = " << current_vals << endl;
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

Eigen::RowVector3d FoldingAnglePositionalConstraintsBuilder::rotate_vec(const Eigen::RowVector3d& pt, const Eigen::RowVector3d& center, const Eigen::Vector3d& axis,
														 double angle) {
	Eigen::Matrix3d rot; rot =  Eigen::AngleAxis<double>(angle, axis);
    return (rot*((pt-center).transpose())).transpose()+center;
}