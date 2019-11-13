#include "FineCoarseConversion.h"

FineCoarseConversion::FineCoarseConversion(const Dog& fine_dog, const Dog& coarse_dog){
	const DogEdgeStitching& fine_es = fine_dog.getEdgeStitching();
	const DogEdgeStitching& coarse_es = coarse_dog.getEdgeStitching();
	const Eigen::MatrixXd& fine_V = fine_dog.getV();
	const Eigen::MatrixXd& coarse_V = coarse_dog.getV();

	int num_curves = fine_es.stitched_curves.size();
	if(num_curves != coarse_es.stitched_curves.size()){
		cout << "Error: Coarse Dog has different number of stitched curves\n";
		return;
	}

	//ftc_curve and coarse_to_fine for curve edge points on both meshes
	ftc_curve.resize(num_curves);
	ctf_curve.resize(num_curves);
	for(int i=0; i<num_curves; ++i){
		ftc_curve[i].resize(fine_es.stitched_curves[i].size());
		ctf_curve[i].resize(coarse_es.stitched_curves[i].size());
		int coarse_j = 0;
		Eigen::RowVector3d coarse_p = coarse_es.stitched_curves[i][coarse_j].getPositionInMesh(coarse_V);
		for(int fine_j=0; fine_j<fine_es.stitched_curves[i].size(); ++fine_j){
			Eigen::RowVector3d fine_p = fine_es.stitched_curves[i][fine_j].getPositionInMesh(fine_V);
			//cout << "ftc: ["<<i<<"]["<<fine_j<<"] is "<< fine_p[0]<<","<< fine_p[1]<<","<< fine_p[2]<<" vs "<< coarse_p[0]<<","<< coarse_p[1]<<","<< coarse_p[2]<<"\n";
			if( (fine_p - coarse_p).squaredNorm() < 1e-5 ){
				//These points have same coordinates, so they're assumed to be the same
				ctf_curve[i][coarse_j] = fine_j;
				ftc_curve[i][fine_j] = coarse_j++;
				if(coarse_j < coarse_es.stitched_curves[i].size())
				  coarse_p = coarse_es.stitched_curves[i][coarse_j].getPositionInMesh(coarse_V);
			} else {
				ftc_curve[i][fine_j] = -1;
			}
//			cout << "ftc: ["<<i<<"]["<<fine_j<<"] to "<<ftc_curve[i][fine_j]<<"\n";
		}
	}


	ctf_curve_vertices.resize(num_curves);
	ctf_curve_weights.resize(num_curves);
	for(int i=0; i<num_curves; ++i){
		ctf_curve_vertices[i] = -1 * Eigen::MatrixXi::Ones(ftc_curve[i].size(), 3);
		ctf_curve_weights[i] = Eigen::MatrixXd::Zero(ftc_curve[i].size(), 3);

		int j_b = 0;//coarse point before
		for(int j=0; j<ftc_curve[i].size(); ++j){
			if(ftc_curve[i][j] == -1){
				//This is an edge point with no coarse equivalent: Let's approximate it
				// in the end it's a linear combination of three coarse vertices (that might be ... bad)
				int vb1,vb2,va1,va2;//v1,v2 , before,after
				vb1 = coarse_es.stitched_curves[i][j_b].edge.v1;
				vb2 = coarse_es.stitched_curves[i][j_b].edge.v2;
				va1 = coarse_es.stitched_curves[i][j_b+1].edge.v1;
				va2 = coarse_es.stitched_curves[i][j_b+1].edge.v2;

				// the target
				Eigen::RowVector3d fine_p = fine_es.stitched_curves[i][j].getPositionInMesh(fine_V);
				// to be expressed as a weighted sum of these
				Eigen::RowVector3d p0 = coarse_V.row(vb1);
				Eigen::RowVector3d p1 = coarse_V.row(vb2);
				Eigen::RowVector3d p2 = coarse_V.row(va1);
				Eigen::RowVector3d p3 = coarse_V.row(va2);

				Eigen::VectorXi parts(4); parts << vb1, vb2, va1, va2;
				Eigen::VectorXd weights(4);

				//Quite possibly two of the vertices involved are the same
				// if so, put one duplicate at the end
				int eq_b = -1, eq_a = -1;
				if(vb1==va1){ eq_b=0; eq_a=2; swap(p2, p3); swap(va1, va2);}
				else if(vb1==va2){ eq_b=0; eq_a=3; }
				else if(vb2==va1){ eq_b=1; eq_a=2; swap(p2, p3); swap(va1, va2);}
				else if(vb2==va2){ eq_b=1; eq_a=3; }

 				Eigen::VectorXd b(3); b << fine_p[0], fine_p[1], 1.0;
				if(eq_b != -1){
					Eigen::MatrixXd A(3,3);
					A << p0[0], p1[0], p2[0],
					     p0[1], p1[1], p2[1],
							   1.0,   1.0,   1.0;
					weights << A.householderQr().solve(b) , 0.0;
				} else {
					Eigen::MatrixXd A(3,4);
					A << p0[0], p1[0], p2[0], p3[0],
					     p0[1], p1[1], p2[1], p3[1],
							   1.0,   1.0,   1.0,   1.0;
					weights = A.householderQr().solve(b);
				}
//				cout << "ftc: vb1,vb2,va1,va2 are "<<vb1<<","<<vb2<<","<<va1<<","<<va2<<"\n";
//				cout << " ftc ctf: ["<<i<<"]["<<j<<"] is "<<weights[0]<<", "<<weights[1]<<", "<<weights[2]<<", "<<weights[3]<<"\n";
				ctf_curve_vertices[i].row(j) << vb1, vb2, va1;
				ctf_curve_weights[i].row(j) << weights[0], weights[1], weights[2];
			} else {
				j_b = ftc_curve[i][j];
				// Just use the corresponding coarse edge point
				EdgePoint coarse_ep = coarse_es.stitched_curves[i][j_b];
				ctf_curve_vertices[i].row(j) << coarse_ep.edge.v1 , coarse_ep.edge.v2, 0;
				ctf_curve_weights[i].row(j) << coarse_ep.t, 1.0-coarse_ep.t, 0;
			}
		}
	}

	//TEST output
	for(int i=0; i<fine_es.stitched_curves[0].size(); ++i){
		Eigen::RowVector3d approximate = get_fine_curve_approx(coarse_V, 0,i);
		Eigen::RowVector3d genuine = fine_es.stitched_curves[0][i].getPositionInMesh(fine_V);
		cout << "comparison: approx "<<approximate <<"\tvs genuine " <<genuine<<"\n";
	}
}

//TODO implement this
int FineCoarseConversion::fine_to_coarse(int fine){
	return -1;
}

//TODO implement this
int FineCoarseConversion::coarse_to_fine(int coarse){
	return -1;
}

//TODO implement this
Eigen::RowVector3d FineCoarseConversion::get_fine_approx(const Eigen::MatrixXd& coarse_V, int fine_v){
	return Eigen::RowVector3d(0,0,0);
}

Eigen::RowVector3d FineCoarseConversion::get_fine_curve_approx(const Eigen::MatrixXd& coarse_V, int curve_idx, int ep_idx){
	Eigen::RowVector3d res;
	res = coarse_V.row(ctf_curve_vertices[curve_idx](ep_idx,0)) * ctf_curve_weights[curve_idx](ep_idx, 0) +
	  coarse_V.row(ctf_curve_vertices[curve_idx](ep_idx, 1)) * ctf_curve_weights[curve_idx](ep_idx, 1) +
		coarse_V.row(ctf_curve_vertices[curve_idx](ep_idx, 2)) * ctf_curve_weights[curve_idx](ep_idx, 2);
	return res;
}
