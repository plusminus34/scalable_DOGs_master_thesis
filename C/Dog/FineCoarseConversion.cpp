#include "FineCoarseConversion.h"

#include "../GeometryPrimitives/Curve.h"

FineCoarseConversion::FineCoarseConversion(const Dog& fine_dog, const Dog& coarse_dog){
	const DogEdgeStitching& fine_es = fine_dog.getEdgeStitching();
	const DogEdgeStitching& coarse_es = coarse_dog.getEdgeStitching();
	const Eigen::MatrixXd& fine_V = fine_dog.getV();
	const Eigen::MatrixXd& coarse_V = coarse_dog.getV();
	const QuadTopology& fine_qt = fine_dog.getQuadTopology();
	const QuadTopology& coarse_qt = coarse_dog.getQuadTopology();

	ftc = Eigen::VectorXi::Constant(fine_V.rows(), -1);
	ctf = Eigen::VectorXi::Constant(coarse_V.rows(), -1);

	//Start by fixing vertex 0 to be the same in both Dogs
	int fine_origin = 0; int coarse_origin = 0;
	ftc(fine_origin) = coarse_origin;
	ctf(coarse_origin) = fine_origin;
	//TODO check that this is okay

	// Need to know about vertices near the patch boundary in the flooding
	int num_patches = fine_dog.get_submesh_n();
	Eigen::MatrixXi fine_vertex_duplicates = Eigen::MatrixXi::Constant(fine_dog.get_v_num(), num_patches, -1);
	Eigen::MatrixXi coarse_vertex_duplicates = Eigen::MatrixXi::Constant(coarse_dog.get_v_num(), num_patches, -1);
	//vertex_duplicates(i,j)=v means that v is i's duplicate in submesh j
	for(int i=0; i<fine_es.edge_const_1.size(); ++i){
		int va1 = fine_es.edge_const_1[i].v1;
		int va2 = fine_es.edge_const_1[i].v2;
		int patch_a = fine_dog.v_to_submesh_idx(va1);
		int vb1 = fine_es.edge_const_2[i].v1;
		int vb2 = fine_es.edge_const_2[i].v2;
		int patch_b = fine_dog.v_to_submesh_idx(vb1);
		fine_vertex_duplicates(va1, patch_b) = vb1;
		fine_vertex_duplicates(va2, patch_b) = vb2;
		fine_vertex_duplicates(vb1, patch_a) = va1;
		fine_vertex_duplicates(vb2, patch_a) = va2;
	}
	for(int i=0; i<coarse_es.edge_const_1.size(); ++i){
		int va1 = coarse_es.edge_const_1[i].v1;
		int va2 = coarse_es.edge_const_1[i].v2;
		int patch_a = coarse_dog.v_to_submesh_idx(va1);
		int vb1 = coarse_es.edge_const_2[i].v1;
		int vb2 = coarse_es.edge_const_2[i].v2;
		int patch_b = coarse_dog.v_to_submesh_idx(vb1);
		coarse_vertex_duplicates(va1, patch_b) = vb1;
		coarse_vertex_duplicates(va2, patch_b) = vb2;
		coarse_vertex_duplicates(vb1, patch_a) = va1;
		coarse_vertex_duplicates(vb2, patch_a) = va2;
	}
	vector<vector<int>> patch_adj = fine_dog.get_submesh_adjacency();
	//Also needed vertex-vertex adjacency
	vector<vector<int>> fine_vv(fine_V.rows());
	vector<vector<int>> coarse_vv(coarse_V.rows());
	for(int i=0; i<fine_vv.size(); ++i) fine_vv[i].clear();
	for(int i=0; i<coarse_vv.size(); ++i) coarse_vv[i].clear();
	for(int i=0; i<fine_qt.E.rows(); ++i){
		fine_vv[fine_qt.E(i,0)].push_back(fine_qt.E(i,1));
		fine_vv[fine_qt.E(i,1)].push_back(fine_qt.E(i,0));
	}
	for(int i=0; i<coarse_qt.E.rows(); ++i){
		coarse_vv[coarse_qt.E(i,0)].push_back(coarse_qt.E(i,1));
		coarse_vv[coarse_qt.E(i,1)].push_back(coarse_qt.E(i,0));
	}

	//There will be a set of vertices to consider, until it finishes
	const int UNCHECKED = -1;
	const int FINEONLY_I = -2;//is between two coarse points
	const int UNDECIDED = -3;
	const int WILLBECHECKED = -4;
	const int FINEONLY_X = -5;//has only fineonly neighbours
	const int COARSEONLY = -6;
	const int LINK = -7;
	//Start by checking the neighbours of the origin
	vector<int> tocheck(0);
	vector<int> tocheck_next = fine_vv[fine_origin];
	vector<int> coarse_tocheck(0);
	vector<int> coarse_tocheck_next = coarse_vv[coarse_origin];
	int flooding_iteration = 0;
	vector<bool> patch_reached(fine_dog.get_submesh_n(), false);
	patch_reached[fine_dog.v_to_submesh_idx(fine_origin)] = true;

	//Flooding to determine type of all vertices
	while(tocheck_next.size() > 0){
		tocheck = tocheck_next;
		tocheck_next.clear();
		++flooding_iteration;
		//cout << "flood: starting iteration "<<flooding_iteration<<"\n";

		if(flooding_iteration%2 == 1){
			//Every odd iteration has no link points
			for(int i=0; i<tocheck.size(); ++i){
				ftc(tocheck[i]) = FINEONLY_I;
				//cout << "flood: "<<tocheck[i]<<" is fine-only and has a coarse neighbour\n";
			}
		} else {
			//Even iterations require more thought
			for(int i=0; i<tocheck.size(); ++i){
				ftc(tocheck[i]) = FINEONLY_X;
				//cout << "flood: "<<tocheck[i]<<" is an x unless stated otherwise\n";
			}
			//Do a coarse iteration
			coarse_tocheck = coarse_tocheck_next;
			coarse_tocheck_next.clear();
			for(int coarse_i=0; coarse_i<coarse_tocheck.size(); ++coarse_i){
				int coarse_u = coarse_tocheck[coarse_i];
				int coarse_patch = coarse_dog.v_to_submesh_idx(coarse_u);
				Eigen::RowVector3d coarse_p = coarse_V.row(coarse_u);
				for(int fine_i=0; fine_i<tocheck.size(); ++fine_i){
					int fine_u = tocheck[fine_i];
					int fine_patch = fine_dog.v_to_submesh_idx(fine_u);
					Eigen::RowVector3d fine_p = fine_V.row(fine_u);
					//Compare the current vertices to be checked (coarse and fine) to find link points
					if(coarse_patch==fine_patch && (coarse_p-fine_p).squaredNorm() < 0.0001){//TODO better check for equality?
						// They are the same
						ctf(coarse_u) = fine_u;
						ftc(fine_u) = coarse_u;
						cout << "flood: linking coarse "<<coarse_u<<" to fine "<<fine_u <<"\n";


						//Try to cross the border between submeshes
						for(int i=0; i<patch_adj[fine_patch].size(); ++i){
							int other_patch = patch_adj[fine_patch][i];
							if( !patch_reached[other_patch] && fine_vertex_duplicates(fine_u, other_patch) != -1){
								// There are a duplicate of vertices fine_u,coarse_u on other_patch
								patch_reached[other_patch] = true;
								int other_fine_u = fine_vertex_duplicates(fine_u, other_patch);
								int other_coarse_u = coarse_vertex_duplicates(coarse_u, other_patch);
								ftc(other_fine_u) = other_coarse_u;
								ctf(other_coarse_u) = other_fine_u;
								// Check neighbours on the other patch in the next iteration
								for(int j=0; j<fine_vv[other_fine_u].size(); ++j){
									int willcheck = fine_vv[other_fine_u][j];
									ftc(willcheck) = WILLBECHECKED;
									tocheck_next.push_back(willcheck);
								}
								for(int j=0; j<coarse_vv[other_coarse_u].size(); ++j){
									int willcheck = coarse_vv[other_coarse_u][j];
									ctf(willcheck) = WILLBECHECKED;
									coarse_tocheck_next.push_back(willcheck);
								}
							}
						}

					}
				}
				for(int j=0; j<coarse_vv[coarse_tocheck[coarse_i]].size(); ++j){
					int coarse_u = coarse_vv[coarse_tocheck[coarse_i]][j];
					if(ctf(coarse_u) == UNCHECKED){
						ctf(coarse_u) = WILLBECHECKED;
						coarse_tocheck_next.push_back(coarse_u);
//						cout << "flood: coarse["<<coarse_u<<"] must be checked next coarse iteration\n";
					}
				}
			}
		}

		//Set unchecked neighbours to be checked in the next iteration
		for(int i=0; i<tocheck.size(); ++i){
			for(int j=0; j<fine_vv[tocheck[i]].size(); ++j){
				int v = fine_vv[tocheck[i]][j];
				if(ftc(v) == UNCHECKED){
					ftc(v) = WILLBECHECKED;
					tocheck_next.push_back(v);
	//				cout << "flood: "<<v<<" must be checked next iteration\n";
				}
			}
		}
	}

	// The rest is curves and some output later
	int num_curves = fine_es.stitched_curves.size();
	if(num_curves != coarse_es.stitched_curves.size()){
		cout << "Error: Coarse Dog has different number of stitched curves\n";
		return;
	}

	//fine_to_coarse and coarse_to_fine for curve edge points on both meshes
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
				cout << "ctf_curve "<<i<<" "<<coarse_j<<" of "<<ctf_curve[i].size()<<"   is now "<<fine_j<<"\n";
				ftc_curve[i][fine_j] = coarse_j++;
				if(coarse_j < coarse_es.stitched_curves[i].size())
				  coarse_p = coarse_es.stitched_curves[i][coarse_j].getPositionInMesh(coarse_V);
			} else {
				ftc_curve[i][fine_j] = -1;
			}
//			cout << "ftc: ["<<i<<"]["<<fine_j<<"] to "<<ftc_curve[i][fine_j]<<"\n";
		}
	}

	// Preparing approximation of fine curve points using coarse ones
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

	// Find curve offsets for interpolating from coarse curve
	vector<SurfaceCurve> coarse_curves(coarse_es.stitched_curves.size());
	for(int i=0; i<coarse_curves.size(); ++i)
		coarse_curves[i].edgePoints = coarse_es.stitched_curves[i];
	ctf_curve_offsets.resize(coarse_curves.size());
	for(int i=0; i<coarse_curves.size(); ++i){
		SurfaceCurve fine_s_curve;
		fine_s_curve.edgePoints = fine_es.stitched_curves[i];
		Curve coarse_curve(coarse_curves[i], coarse_V);
		Curve fine_curve(fine_s_curve, fine_V);
		ctf_curve_offsets[i].resize(coarse_curve.len.size());
		int current_coarse = 0;
		int next_coarse = 1;
		ctf_curve_offsets[i][0].push_back(0.0);
		while(next_coarse < coarse_curves[i].edgePoints.size()){
			int fine_begin = ctf_curve[i][current_coarse];
			int fine_end = ctf_curve[i][next_coarse];
			if(fine_begin < fine_end){
				int n_steps = fine_end - fine_begin;
				vector<double> sum_len(n_steps);
				sum_len[0] = fine_curve.len[fine_begin];
				for(int j=1; j < n_steps; ++j){
					sum_len[j] = sum_len[j-1] + fine_curve.len[j+fine_begin];
				}
				double total_fine_len = sum_len[sum_len.size()-1];
				for(int j=0; j < n_steps; ++j){
					double offset = sum_len[j] / total_fine_len;
					ctf_curve_offsets[i][current_coarse].push_back(offset);
				}
			}
			current_coarse = next_coarse;
			++next_coarse;
		}
	}
	      /*
	      for(int i=0;i<offsets.size(); ++i)
	        for(int j=0;j<offsets[i].size(); ++j)
	          for(int k=0;k<offsets[i][j].size(); ++k)
	            cout <<"offsets["<<i<<"]["<<j<<"]["<<k<<"] =\t"<<offsets[i][j][k]<<endl;
	            */

	//TEST output
/*
	for(int i=0; i<fine_es.stitched_curves[0].size(); ++i){
		Eigen::RowVector3d approximate = get_fine_curve_approx(coarse_V, 0,i);
		Eigen::RowVector3d genuine = fine_es.stitched_curves[0][i].getPositionInMesh(fine_V);
		cout << "comparison: approx "<<approximate <<"\tvs genuine " <<genuine<<"\n";
	}*/
	print();

}

void FineCoarseConversion::print(){
	int links = 0;
	int fineonly = 0;
	int coarseonly = 0;
	for(int i=0; i<ftc.size(); ++i){
		cout << "fine_to_coarse: "<<i;
		if(ftc[i] > -1) {
			++links;
			cout <<" links to " << ftc[i] << " and accordingly coarse_to_fine "<<ftc[i]<<" is "<<ctf[ftc[i]]<<"\n";
		}else{
			++fineonly;
			cout <<" is "<<ftc[i]<<"\n";
		}
	}
	for(int i=0; i<ctf.size(); ++i){
		if(ctf[i]<0){
			++coarseonly;
			cout << "coarse_to_fine: coarse "<<i<<" is sad and alone "<<ctf[i]<<"\n";
		}else{
			cout << "coarse_to_fine: coarse "<<i<<" links to "<<ctf[i]<<" and accordingly fine_to_coarse "<<ctf[i]<<" is "<<ftc[ctf[i]]<<"\n";
		}
	}
	cout << "From "<<ftc.size()<<" fine vertices and "<<ctf.size()<<" coarse vertices:\n";
	cout << "There are "<<links<<" link vertices, "<<fineonly<<" fine-only and "<<coarseonly<<" coarse-only vertices.\n";
}

int FineCoarseConversion::fine_to_coarse(int fine){
	return ftc(fine);
}

int FineCoarseConversion::coarse_to_fine(int coarse){
	return ctf(coarse);
}

int FineCoarseConversion::coarse_to_fine_curve(int curve_idx, int ep_idx){
	return ctf_curve[curve_idx][ep_idx];
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
