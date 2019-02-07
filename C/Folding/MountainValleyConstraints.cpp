#include "MountainValleyConstraints.h"

using namespace std;

MountainValleyConstraints::MountainValleyConstraints(const Dog& dog,const DogEdgeStitching& edgeStitching, 
																				  std::vector<bool> is_mountain) :
																					init_dog(dog), eS(edgeStitching), is_mountain(is_mountain) {
    // TODO: different handling of vertices (maybe no consatraints there?)
	const_n= 2*edgeStitching.edge_coordinates.size();
	IJV.resize(15*const_n);

	// initialize flip_binormal
	flip_binormal.resize(eS.stitched_curves.size());
	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];	
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			// get the 3 points
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			auto ep_p = ep.getPositionInMesh(init_dog.getV()); auto ep_b_p = ep_b.getPositionInMesh(init_dog.getV()); auto ep_f_p = ep_f.getPositionInMesh(init_dog.getV());

			Eigen::RowVector3d b = (ep_f_p-ep_p).cross(ep_p-ep_b_p).normalized();
			Eigen::RowVector3d z_axis = Eigen::RowVector3d(0,0,1);
			// if closer to -1,0,0 than to 1,0,0 mark flip_binormals as "yes"
			// necessary to define M/V assignment correctly
			flip_binormal[curve_i].push_back( (b.dot(z_axis) < 0 ) );
		}
	}
}

Eigen::VectorXd MountainValleyConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
  	
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			bool flip_binormal_sign = (is_mountain[curve_i] ^ flip_binormal[curve_i][edge_idx-1]);
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if (flip_binormal_sign) {
				// Flip the edge direction
				std::swap(ep_b,ep_f);
			}
			int fold_v_indices[2]; init_dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);

			int ep_0_v1_i(ep.edge.v1), ep_0_v2_i(ep.edge.v2); const double ep_0_t(ep.t);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
			const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			for (int const_side_i = 0; const_side_i < 2 ; const_side_i++) {
				int folded_v_i = fold_v_indices[const_side_i];
				const double folded_v_x(x(folded_v_i)); const double folded_v_y(x(folded_v_i+1*vnum)); const double folded_v_z(x(folded_v_i+2*vnum));

				double t2 = ep_0_t-1.0;
				double t3 = ep_0_t*ep_0_v1_y;
				double t4 = ep_b_t-1.0;
				double t5 = ep_0_t*ep_0_v1_x;
				double t6 = ep_f_t-1.0;
				double t7 = ep_b_v2_x*t4;
				double t10 = ep_0_v2_x*t2;
				double t8 = t5+t7-t10-ep_b_t*ep_b_v1_x;
				double t9 = ep_0_v2_z*t2;
				double t11 = ep_f_v2_x*t6;
				double t12 = ep_0_v2_y*t2;
				double t13 = ep_b_v2_y*t4;
				double t14 = ep_f_t*ep_f_v1_z;
				double t16 = ep_0_t*ep_0_v1_z;
				double t15 = t9+t14-t16-ep_f_v2_z*t6;
				double t17 = ep_b_t*ep_b_v1_z;
				double t18 = ep_f_v2_y*t6;
				constVals(const_cnt++) = -(t8*(t3+t18-ep_f_t*ep_f_v1_y-ep_0_v2_y*t2)-(t3+t13-ep_b_t*ep_b_v1_y-ep_0_v2_y*t2)*(t5+t11-ep_f_t*ep_f_v1_x-ep_0_v2_x*t2))*(folded_v_z+t9-ep_0_t*ep_0_v1_z)-(t8*t15-(t5-t10+t11-ep_f_t*ep_f_v1_x)*(t9+t17-ep_0_t*ep_0_v1_z-ep_b_v2_z*t4))*(folded_v_y-t3+t12)+(t15*(t3-t12+t13-ep_b_t*ep_b_v1_y)-(t3-t12+t18-ep_f_t*ep_f_v1_y)*(t9-t16+t17-ep_b_v2_z*t4))*(folded_v_x-t5+t10);
			}
		}
	}
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  return constVals;
}


void MountainValleyConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {

	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0; int ijv_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			bool flip_binormal_sign = (is_mountain[curve_i] ^ flip_binormal[curve_i][edge_idx-1]);
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if (flip_binormal_sign) {
				// Flip the edge direction
				std::swap(ep_b,ep_f);
			}
			int fold_v_indices[2]; init_dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);

			int ep_0_v1_i(ep.edge.v1), ep_0_v2_i(ep.edge.v2); const double ep_0_t(ep.t);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
			const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			for (int const_side_i = 0; const_side_i < 2 ; const_side_i++) {
				int folded_v_i = fold_v_indices[const_side_i];
				const double folded_v_x(x(folded_v_i)); const double folded_v_y(x(folded_v_i+1*vnum)); const double folded_v_z(x(folded_v_i+2*vnum));

				double t2 = ep_0_t-1.0;
				double t3 = ep_0_t*ep_0_v1_z;
				double t4 = ep_b_t-1.0;
				double t5 = ep_0_t*ep_0_v1_y;
				double t6 = ep_f_t-1.0;
				double t7 = ep_b_v2_y*t4;
				double t9 = ep_0_v2_y*t2;
				double t28 = ep_b_t*ep_b_v1_y;
				double t8 = t5+t7-t9-t28;
				double t10 = ep_f_v2_y*t6;
				double t11 = ep_0_v2_z*t2;
				double t12 = ep_b_v2_z*t4;
				double t13 = ep_f_v2_z*t6;
				double t23 = ep_f_t*ep_f_v1_z;
				double t14 = t3-t11+t13-t23;
				double t22 = ep_b_t*ep_b_v1_z;
				double t15 = t3-t11+t12-t22;
				double t16 = ep_0_t*ep_0_v1_x;
				double t17 = ep_b_v2_x*t4;
				double t19 = ep_0_v2_x*t2;
				double t26 = ep_b_t*ep_b_v1_x;
				double t18 = t16+t17-t19-t26;
				double t20 = ep_f_v2_x*t6;
				double t21 = folded_v_z-t3+t11;
				double t24 = ep_0_t*t14;
				double t25 = t24-ep_0_t*t15;
				double t35 = ep_f_t*ep_f_v1_y;
				double t27 = t5-t9+t10-t35;
				double t31 = ep_f_t*ep_f_v1_x;
				double t29 = t16-t19+t20-t31;
				double t30 = ep_0_t*t18;
				double t32 = t30-ep_0_t*t29;
				double t33 = folded_v_y-t5+t9;
				double t34 = ep_0_t*t8;
				double t36 = t34-ep_0_t*t27;
				double t37 = folded_v_x-t16+t19;
				double t38 = t2*t14;
				double t39 = t14*t18;
				double t40 = t2*t18;
				double t41 = t40-t2*t29;
				double t42 = t2*t8;
				double t43 = t42-t2*t27;
				double t44 = t8*t29;
				double t45 = t8*t14;
				double t46 = t39-t15*t29;
				double t47 = t44-t18*t27;

				//[ ep_0_v1_x, ep_0_v1_y, ep_0_v1_z, ep_0_v2_x, ep_0_v2_y, ep_0_v2_z, ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, folded_v_x, folded_v_y, folded_v_z]
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v1_i,t21*t36+t25*t33+ep_0_t*(t8*(t3+t13-ep_f_t*ep_f_v1_z-ep_0_v2_z*t2)-(t3+t12-ep_b_t*ep_b_v1_z-ep_0_v2_z*t2)*(t5+t10-ep_f_t*ep_f_v1_y-ep_0_v2_y*t2)));
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v1_i+vnum,-ep_0_t*(t39-t15*(t16+t20-ep_f_t*ep_f_v1_x-ep_0_v2_x*t2))-t21*t32-t25*t37);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v1_i+2*vnum,-ep_0_t*t47+t32*t33-t36*t37);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v2_i,-t21*t43-t33*(t38-t2*t15)-t2*(t45-t15*t27));
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v2_i+vnum,t2*t46+t21*t41+t37*(t38-t2*t15));
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_0_v2_i+2*vnum,t2*t47-t33*t41+t37*t43);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i,-ep_b_t*t14*t33+ep_b_t*t21*t27);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+vnum,-ep_b_t*t21*t29+ep_b_t*t14*t37);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+2*vnum,ep_b_t*t29*t33-ep_b_t*t27*t37);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i,t4*t14*t33-t4*t21*t27);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+vnum,t4*t21*t29-t4*t14*t37);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+2*vnum,-t4*t29*t33+t4*t27*t37);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i,-ep_f_t*t8*t21+ep_f_t*t15*t33);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+vnum,ep_f_t*t18*t21-ep_f_t*t15*t37);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+2*vnum,ep_f_t*t8*t37-ep_f_t*t18*t33);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i,t6*t8*t21-t6*t15*t33);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+vnum,-t6*t18*t21+t6*t15*t37);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+2*vnum,-t6*t8*t37+t6*t18*t33);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,folded_v_i,-t45+t15*t27);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,folded_v_i+vnum,t46);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,folded_v_i+2*vnum,t47);

				const_cnt++;
			}
		}
	}
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}