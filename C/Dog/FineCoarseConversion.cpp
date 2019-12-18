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
	ftc_edge = Eigen::VectorXi::Constant(fine_V.rows(), -1);

	//Usually vertex 0 is the same in both Dogs
	int fine_origin = 0; int coarse_origin = 0;
	// sometimes it isn't
	while( (fine_V.row(fine_origin) - coarse_V.row(coarse_origin) ).squaredNorm() > 0.0001
			|| fine_dog.v_to_submesh_idx(fine_origin) != coarse_dog.v_to_submesh_idx(coarse_origin)){
		++fine_origin;
		if(fine_origin == fine_V.rows()){
			fine_origin = 0;
			++coarse_origin;
			if(coarse_origin == coarse_V.rows()){
				coarse_origin = 0;
				cout << "Error: Couldn't find suitable origin for fine-to-coarse conversion\n";
				break;
			}
		}
	}
	ftc(fine_origin) = coarse_origin;
	ctf(coarse_origin) = fine_origin;

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
	//Also needed: vertex-vertex adjacency
	vector<vector<int>> fine_vv(fine_V.rows());
	vector<vector<int>> coarse_vv(coarse_V.rows());
	map<pair<int,int>, int> fine_edge_to_row;
	map<pair<int,int>, int> coarse_edge_to_row;
	for(int i=0; i<fine_vv.size(); ++i) fine_vv[i].clear();
	for(int i=0; i<coarse_vv.size(); ++i) coarse_vv[i].clear();
	for(int i=0; i<fine_qt.E.rows(); ++i){
		int v1 = fine_qt.E(i,0);
		int v2 = fine_qt.E(i,1);
		fine_vv[v1].push_back(v2);
		fine_vv[v2].push_back(v1);
		fine_edge_to_row[pair<int,int>(v1,v2)] = i;
	}
	for(int i=0; i<coarse_qt.E.rows(); ++i){
		int v1 = coarse_qt.E(i,0);
		int v2 = coarse_qt.E(i,1);
		coarse_vv[v1].push_back(v2);
		coarse_vv[v2].push_back(v1);
		coarse_edge_to_row[pair<int,int>(v1,v2)] = i;
	}

	//There will be a set of vertices to consider, until it finishes
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
			/*
			TODO find out why exactly coarse flooding fails (brute force method used until then)
			//Do a coarse iteration
			coarse_tocheck = coarse_tocheck_next;
			coarse_tocheck_next.clear();
			for(int coarse_i=0; coarse_i<coarse_tocheck.size(); ++coarse_i){
				int coarse_u = coarse_tocheck[coarse_i];
				*/
			for(int coarse_i = 0; coarse_i<coarse_V.rows(); ++coarse_i){
				int coarse_u = coarse_i;
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
//									coarse_tocheck_next.push_back(willcheck);
								}
							}
						}

					}
				}
				/*
				for(int j=0; j<coarse_vv[coarse_tocheck[coarse_i]].size(); ++j){
					int coarse_u = coarse_vv[coarse_tocheck[coarse_i]][j];
					if(ctf(coarse_u) == UNCHECKED){
						ctf(coarse_u) = WILLBECHECKED;
						coarse_tocheck_next.push_back(coarse_u);
//						cout << "flood: coarse["<<coarse_u<<"] must be checked next coarse iteration\n";
					}
				}
				*/
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
				//cout << "ctf_curve "<<i<<" "<<coarse_j<<" of "<<ctf_curve[i].size()<<"   is now "<<fine_j<<"\n";
				ftc_curve[i][fine_j] = coarse_j++;
				if(coarse_j < coarse_es.stitched_curves[i].size())
				  coarse_p = coarse_es.stitched_curves[i][coarse_j].getPositionInMesh(coarse_V);
			} else {
				ftc_curve[i][fine_j] = -1;
			}
//			cout << "ftc: ["<<i<<"]["<<fine_j<<"] to "<<ftc_curve[i][fine_j]<<"\n";
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


	// Here comes fine_to_coarse_edge plus data to update COARSEONLY later
	coarseonly_adjacent_links.resize(ctf.size());
	coarseonly_adjacent_fineonly.resize(ctf.size());
	for(int i=0; i<ctf.size(); ++i){
		coarseonly_adjacent_links[i].clear();
		coarseonly_adjacent_fineonly[i].clear();
	}
	//get list of FINEONLY_I vertices
	vector<int> I_vertex(0);
	for(int i=0; i<ftc.size(); ++i) if(ftc(i) == FINEONLY_I) I_vertex.push_back(i);

	for(int i=0; i<I_vertex.size(); ++i){
		int u = I_vertex[i];
		int first_link_nb = -1;
		for(int j=0; j<fine_vv[u].size(); ++j){
			int v = fine_vv[u][j];
			if(ftc(v) > -1) {
				if(first_link_nb == -1) first_link_nb = v;
				else {
					// found second LINK neighbour
					int coarse_nb_1 = ftc(first_link_nb);
					int coarse_nb_2 = ftc(v);
					if(coarse_nb_1 > coarse_nb_2) swap(coarse_nb_1, coarse_nb_2);
					ftc_edge(u) = coarse_edge_to_row.at(pair<int,int>(coarse_nb_1, coarse_nb_2));
					break;
				}
			}
		}
		// there is no second LINK neighbour
		if(ftc_edge(u)>-1) continue;
		int coarse_nb_1 = ftc(first_link_nb);
		Eigen::RowVector3d p = coarse_V.row(coarse_nb_1);
		Eigen::RowVector3d a = fine_V.row(u);
		// have to check for a suited COARSEONLY neighbour
		for(int k=0; k<coarse_vv[coarse_nb_1].size(); ++k){
			int coarse_v = coarse_vv[coarse_nb_1][k];
			if(ctf(coarse_v) > -1) continue;
			Eigen::RowVector3d b = coarse_V.row(coarse_v);
			double dot = (a-p).dot(b-p);
			if (dot>0){
				//add data for COARSEONLY update
				coarseonly_adjacent_fineonly[coarse_v].push_back(u);
				coarseonly_adjacent_links[coarse_v].push_back(first_link_nb);
				//and the edge mapping as expected
				int coarse_nb_2 = coarse_v;
				if(coarse_nb_1 > coarse_nb_2) swap(coarse_nb_1, coarse_nb_2);
				ftc_edge(u) = coarse_edge_to_row.at(pair<int,int>(coarse_nb_1, coarse_nb_2));
				break;
			}
		}
	}

	for(int i=0;i<ctf_curve_offsets.size(); ++i)
		for(int j=0;j<ctf_curve_offsets[i].size(); ++j)
			for(int k=0;k<ctf_curve_offsets[i][j].size(); ++k)
				cout <<"offsets["<<i<<"]["<<j<<"]["<<k<<"] =\t"<<ctf_curve_offsets[i][j][k]<<endl;
			getInterpolatedCurveCoords(fine_dog,coarse_dog,0);
}

Dog FineCoarseConversion::init_from_fine_dog(const Dog& fine_dog){
	//TODO do this
	const int fine_v_num = fine_dog.get_v_num();
	const DogEdgeStitching& fine_es = fine_dog.getEdgeStitching();
	const Eigen::MatrixXd& fine_V = fine_dog.getV();
	const Eigen::MatrixXi& fine_F = fine_dog.getF();
	const QuadTopology& fine_qt = fine_dog.getQuadTopology();
	const int num_submeshes = fine_dog.get_submesh_n();
	const int num_curves = fine_es.stitched_curves.size();

	//very helpful: vertex-vertex adjacency
	vector<vector<int>> fine_vv(fine_v_num);
	for(int i=0; i<fine_v_num; ++i) fine_vv[i].clear();
	for(int i=0; i<fine_qt.E.rows(); ++i){
		int v1 = fine_qt.E(i,0); int v2 = fine_qt.E(i,1);
		fine_vv[v1].push_back(v2); fine_vv[v2].push_back(v1);
	}
	//probably also helpful: crease vertices to duplicates
	vector<vector<int>> fine_duplicates(fine_v_num);
	for(int i=0; i<fine_v_num; ++i) fine_duplicates[i].clear();
	for(int i=0; i<fine_es.edge_const_1.size(); ++i){
		int v1 = fine_es.edge_const_1[i].v1;
		int v2 = fine_es.edge_const_1[i].v2;
		int w1 = fine_es.edge_const_2[i].v1;
		int w2 = fine_es.edge_const_2[i].v2;
		fine_duplicates[v1].push_back(w1);
		fine_duplicates[v2].push_back(w2);
		fine_duplicates[w1].push_back(v1);
		fine_duplicates[w2].push_back(v2);
	}

	// to be built
	ftc.resize(fine_v_num);
	vector<int> fine_label(fine_v_num, UNCHECKED);// might be useful to save this
	//ctf not yet, number of coarse vertices is unknown
	ftc_edge.resize(fine_v_num); ftc_edge.setConstant(-1);
	ftc_curve.resize(num_curves);
	ctf_curve.resize(num_curves);
	entire_coarse_curve_i.resize(num_curves);
	entire_coarse_curve_v.resize(num_curves);
	entire_coarse_curve_w.resize(num_curves);
	for(int i=0; i<num_curves; ++i){
		ftc_curve[i] = vector<int>(fine_es.stitched_curves[i].size(), -1);
		ctf_curve[i].clear();
		entire_coarse_curve_i.clear();
		entire_coarse_curve_v[i].resize(fine_es.stitched_curves[i].size(), 4);
		entire_coarse_curve_w[i].resize(fine_es.stitched_curves[i].size(), 4);
		entire_coarse_curve_v[i].setZero(); entire_coarse_curve_w[i].setZero();
	}
	// ctf_curve_offsets come later
	DogEdgeStitching coarse_es;
	vector<int> coarse_submeshVSize(num_submeshes, 0);
	vector<int> coarse_submeshFSize(num_submeshes, 0);
	vector<vector<int>> submesh_adjacency = fine_dog.get_submesh_adjacency();

	//use something as origin
	int fine_origin = 0;
	for(int i=0;i<fine_v_num;++i){
		if(fine_vv[i].size()==2&&fine_duplicates[i].size()==0){fine_origin=i;break;}
	}
	cout <<"origin: "<<fine_origin<<endl;
	//Spread out from the origin
	// Each fine quad has a LINK vertex, two FINEONLY_I and one FINEONLY_X
	fine_label[fine_origin] = LINK;
	list<int> tocheck;tocheck.push_back(fine_origin);
	while(tocheck.size() > 0){
		int v = tocheck.back();
		tocheck.pop_back();
		vector<int> close_pts(0);
		vector<int> quads = fine_qt.VF[v];
		for(int i=0; i<quads.size(); ++i){
			for(int j=0; j<4; ++j) close_pts.push_back(fine_F(quads[i], j));
		}
		for(int i=0; i<close_pts.size(); ++i){
			int u = close_pts[i];
			if(fine_label[u] != UNCHECKED) continue;
			for(int j=0; j<fine_vv[v].size(); ++j) if(fine_vv[v][j]==u){fine_label[u] = FINEONLY_I; break;}
			if(fine_label[u] != UNCHECKED) continue;
			if(fine_label[v] == LINK) fine_label[u] = FINEONLY_X;
			else if(fine_label[v] == FINEONLY_X) {
				fine_label[u] = LINK;
				// try to cross the patch border
				for(int j=0; j<fine_duplicates[u].size(); ++j){
					int w = fine_duplicates[u][j];
					if(fine_label[w] == UNCHECKED){
						fine_label[w] = LINK;
						tocheck.push_back(w);
					}
				}
			}
			tocheck.push_back(u);
		}
	}
	int links=0, is=0, xs=0, unchecked=0;
	for(int i=0; i<fine_v_num; ++i){
		if(fine_label[i] == FINEONLY_I) ++is;
		else if(fine_label[i] == LINK) ++links;
		else if(fine_label[i] == FINEONLY_X) ++xs;
		else if(fine_label[i] == UNCHECKED) ++unchecked;
	}
	cout << "There are "<<links<<" link vertices, "<<is<<" fineonly_i, "<<xs<<" fineonly_x, "<<unchecked<<" unchecked.\n";
	if(unchecked > 0) cout << "Having unchecked vertices means something is wrong.\n";
	// Labels done, now use them to construct the rest
	// Build coarse V and F
	vector<int> coarse_quad_centers(0);
	ctf.resize(links);//no coarse-only vertices this time!
	Eigen::MatrixXd coarse_V(links, 3);
	int coarse_i=0;
	for(int i=0; i<fine_v_num; ++i){
		if(fine_label[i] == LINK){
			ftc(i) = coarse_i;
			ctf(coarse_i) = i;
			coarse_V.row(coarse_i++) = fine_V.row(i);
			++coarse_submeshVSize[fine_dog.v_to_submesh_idx(i)];
		} else if (fine_label[i] == FINEONLY_X){
			if(fine_qt.VF[i].size() == 4) {
				coarse_quad_centers.push_back(i);
				++coarse_submeshFSize[fine_dog.v_to_submesh_idx(i)];
			}
		}
	}
	cout << "Built coarse_V\n";
	Eigen::MatrixXi coarse_F(coarse_quad_centers.size(), 4);
	for(int i=0; i<coarse_quad_centers.size(); ++i){
		int v = coarse_quad_centers[i];
		for(int j=0; j<4; ++j){
			int quad = fine_qt.VF[v][j];
		  for(int k=0; k<4; ++k){
				int u = fine_F(quad, k);
				if(fine_label[u] == LINK) {coarse_F(i,j) = ftc(u); break;}
			}
		}
	}
	cout << "Built coarse_F\n";
	// Build ftc_edge (mapping FINEONLY_I vertices to a coarse edge)
	QuadTopology coarse_qt; quad_topology(coarse_V, coarse_F, coarse_qt);
	cout << "Built coarse_qt\n";
	for(int i=0; i<coarse_qt.E.rows(); ++i){
		int v1 = ctf(coarse_qt.E(i,0));
		int v2 = ctf(coarse_qt.E(i,1));
		for(int j=0; j<fine_v_num; ++j){
			if(fine_label[j] != FINEONLY_I) continue;
			bool adj1 = false, adj2 = false;
			// find out if v1 and v2 are neighbours of j
			for(int k=0; k<fine_vv[j].size(); ++k){
				if(fine_vv[j][k] == v1) adj1 = true;
				if(fine_vv[j][k] == v2) adj2 = true;
			}
			if(adj1 && adj2){
				ftc_edge(j) = i;
				continue;
			}
		}
	}
	cout << "Built ftc_edge\n";
	// start building coarse_es
	coarse_es.stitched_curves.resize(num_curves);
	for(int i=0; i<num_curves; ++i) coarse_es.stitched_curves[i].clear();
	for(int i=0; i<num_curves; ++i){
		for(int j=0; j<fine_es.stitched_curves[i].size(); ++j){
			int v1,v2,w1,w2;
			EdgePoint ep = fine_es.stitched_curves[i][j];
			fine_dog.get_2_submeshes_vertices_from_edge(ep.edge, v1, v2, w1, w2);
			int c_v1 = -1, c_v2 = -1, c_w1 = -1, c_w2 = -1;
			double c_t = -1.0;
			if(fine_label[v1] == FINEONLY_I && fine_label[v2] == LINK && ftc_edge[v1] > -1){
				c_v1 = coarse_qt.E(ftc_edge[v1], 0);
				c_v2 = coarse_qt.E(ftc_edge[v1], 1);
				c_t = 0.5 * ep.t;
			} else if(fine_label[v2] == FINEONLY_I && fine_label[v1] == LINK && ftc_edge[v2] > -1){
				c_v1 = coarse_qt.E(ftc_edge[v2], 0);
				c_v2 = coarse_qt.E(ftc_edge[v2], 1);
				c_t = 0.5 * ep.t + 0.5;
			}
			if(fine_label[w1] == FINEONLY_I && fine_label[w2] == LINK && ftc_edge[w1] > -1){
				c_w1 = coarse_qt.E(ftc_edge[w1], 0);
				c_w2 = coarse_qt.E(ftc_edge[w1], 1);
				c_t = 0.5 * ep.t;
			} else if(fine_label[w2] == FINEONLY_I && fine_label[w1] == LINK && ftc_edge[w2] > -1){
				c_w1 = coarse_qt.E(ftc_edge[w2], 0);
				c_w2 = coarse_qt.E(ftc_edge[w2], 1);
				c_t = 0.5 * ep.t + 0.5;
			}
			if(c_t < 0) continue;// skip if no coarse submesh contains ep
			//if( (j>0&&j<fine_es.stitched_curves[i].size()-1) && (c_v1 == -1 || c_w1 == -1) ) {continue;}//only take endpoints as not-so-good curve points
			entire_coarse_curve_i[i].push_back(j);
			if(c_v1 > -1){
				entire_coarse_curve_v[i](j, 0) = c_v1;
				entire_coarse_curve_v[i](j, 1) = c_v2;
				entire_coarse_curve_w[i](j, 0) = c_t;
				entire_coarse_curve_w[i](j, 1) = 1.0 - c_t;
			}
			if(c_w1 > -1){
				entire_coarse_curve_v[i](j, 2) = c_w1;
				entire_coarse_curve_v[i](j, 3) = c_w2;
				entire_coarse_curve_w[i](j, 2) = c_t;
				entire_coarse_curve_w[i](j, 3) = 1.0 - c_t;
			}
			ftc_curve[i][j] = COARSEONLY;
			cout << "ftc curve: extra not-so-good coarse point\n";
			if(c_v1 == -1 || c_w1 == -1) continue;// skip unless both coarse submeshes contain ep
			entire_coarse_curve_w[i].row(j) *= 0.5;
			int coarse_j = coarse_es.stitched_curves[i].size();
			coarse_es.stitched_curves[i].push_back( EdgePoint(Edge(c_v1,c_v2), c_t) );
			ftc_curve[i][j] = coarse_j;
			ctf_curve[i].push_back(j);
			cout << "ftc curve: fine "<<i<<" "<<j<<" to coarse "<<i<<" "<<coarse_j<<endl;
		}
	}
	//return Dog(coarse_V, coarse_F);
	coarse_es.submesh_to_edge_pt.resize(num_submeshes);
	coarse_es.submesh_to_bnd_edge.resize(num_submeshes);
	for(int i=0; i<num_submeshes; ++i){
		coarse_es.submesh_to_edge_pt.clear();
		coarse_es.submesh_to_bnd_edge.clear();
	}
	for(int i=0; i<fine_es.edge_const_1.size(); ++i){
		int v1 = fine_es.edge_const_1[i].v1; int v2 = fine_es.edge_const_1[i].v2;
		int w1 = fine_es.edge_const_2[i].v1; int w2 = fine_es.edge_const_2[i].v2;
		int c_ev, c_ew = -1;
		double c_t;
		if(ftc_edge[v1]>-1 && ftc_edge[w1]>-1 && fine_label[v2] == LINK && fine_label[w2] == LINK){
			c_ev = ftc_edge[v1]; c_ew = ftc_edge[w1];
			c_t = 0.5 * fine_es.edge_coordinates[i];
		} else if (ftc_edge[v2]>-1 && ftc_edge[w2]>-1 && fine_label[v1] == LINK && fine_label[w1] == LINK){
			c_ev = ftc_edge[v2]; c_ew = ftc_edge[w2];
			c_t=0.5 * fine_es.edge_coordinates[i] + 0.5;
		}
		if(c_ew < 0) continue;
		int c_v1 = coarse_qt.E(c_ev,0); int c_v2 = coarse_qt.E(c_ev,1);
		int c_w1 = coarse_qt.E(c_ew,0); int c_w2 = coarse_qt.E(c_ew,1);
		Edge edge_v(c_v1, c_v2); Edge edge_w(c_w1, c_w2);
		int idx = coarse_es.edge_const_1.size();

		coarse_es.edge_const_1.push_back(edge_v);
		coarse_es.edge_const_2.push_back(edge_w);
		coarse_es.edge_coordinates.push_back(c_t);

		coarse_es.submesh_to_edge_pt[fine_dog.v_to_submesh_idx(v1)].push_back(idx);
		coarse_es.submesh_to_edge_pt[fine_dog.v_to_submesh_idx(w1)].push_back(idx);
		//These down here assume that no coarse crease point appears 3 or more times
		coarse_es.edge_to_duplicates[edge_v] = idx;
		coarse_es.edge_to_duplicates[edge_w] = idx;
		coarse_es.multiplied_edges_start.push_back(idx);
		coarse_es.multiplied_edges_num.push_back(1);
	}
	cout << "Built edge_const and edge_coordinates for coarse edgeStitching\n";
	cout << "   and edge_to_duplicates, multiplied_edges_start, multiplied_edges_num as well (probably)\n";
	cout << "   actually, the coarse edgeStitching should be ready now\n";
	// Build ctf_curve_offsets
	ctf_curve_offsets.resize(num_curves);
	for(int i=0; i<num_curves; ++i){
		//init offsets[i]
		ctf_curve_offsets[i].resize(entire_coarse_curve_i[i].size()-1);
		for(int j=0; j<ctf_curve_offsets[i].size();++j)ctf_curve_offsets[i][j].clear();
		ctf_curve_offsets[i][0].push_back(0.0);
		//Get coarse and fine curve coordinates
		Eigen::MatrixXd coarse_coords(entire_coarse_curve_i[i].size(), 3);
		coarse_coords.setZero();
		for(int j=0; j<entire_coarse_curve_i[i].size(); ++j){
			int row = entire_coarse_curve_i[i][j];
			for(int k=0; k<4; ++k)
				if(entire_coarse_curve_w[i](row,k)>0)
					coarse_coords.row(j) += entire_coarse_curve_w[i](row,k) * coarse_V.row(entire_coarse_curve_v[i](row,k));
			}
		Eigen::MatrixXd fine_coords(fine_es.stitched_curves[i].size(), 3);
		for(int j=0; j<fine_coords.rows(); ++j)
			fine_coords.row(j) = fine_es.stitched_curves[i][j].getPositionInMesh(fine_V);
		//actual offset computation
		/*
		for(int coarse_j=0; coarse_j < entire_coarse_curve_i[i].size()-1; ++coarse_j){
			int fine_begin = entire_coarse_curve_i[i][coarse_j];
			int fine_end = entire_coarse_curve_i[i][coarse_j+1];
			*/
		for(int coarse_j=0; coarse_j < coarse_es.stitched_curves[i].size()-1; ++coarse_j){
			int fine_begin = ctf_curve[i][coarse_j];
			int fine_end = ctf_curve[i][coarse_j+1];
			if(fine_begin < fine_end){
				int n_steps = fine_end - fine_begin;
				vector<double> sum_len(n_steps);
				sum_len[0] = (fine_coords.row(fine_begin+1) - fine_coords.row(fine_begin)).norm();
				for(int j=1; j < n_steps; ++j){
					sum_len[j] = sum_len[j-1] + (fine_coords.row(fine_begin+j+1) - fine_coords.row(fine_begin+j)).norm();
				}
				double total_fine_len = sum_len[sum_len.size()-1];
				for(int j=0; j < n_steps; ++j){
					double offset = sum_len[j] / total_fine_len;
					ctf_curve_offsets[i][coarse_j].push_back(offset);
				}
			}
		}
	}
	for(int i=0;i<ctf_curve_offsets.size(); ++i)
		for(int j=0;j<ctf_curve_offsets[i].size(); ++j)
			for(int k=0;k<ctf_curve_offsets[i][j].size(); ++k)
				cout <<"offsets["<<i<<"]["<<j<<"]["<<k<<"] =\t"<<ctf_curve_offsets[i][j][k]<<endl;
	cout << "Built offsets for curve interpolation\n";

	// Final preparations
	for(int j=0;j<ftc_curve[0].size(); ++j){
		cout << "ftc curve [0]["<<j<<"] = "<<ftc_curve[0][j]<<"\n";
		EdgePoint ep= fine_es.stitched_curves[0][j];
		cout << "ftc curve and its edgp is "<<ep.edge.v1<<" - "<<ep.edge.v2 <<"    t="<<ep.t<<endl;
	}

	//Eigen::MatrixXd coarse_V_ren = coarse_V;
	Eigen::MatrixXi coarse_F_ren = Fsqr_to_F(coarse_F);
	coarse_es.edge_coordinates_precise.clear();
	for(int i=0; i<ftc.size(); ++i){
		if (fine_label[i]== FINEONLY_I) ftc(i)=FINEONLY_I;
		else if (fine_label[i]== FINEONLY_X) ftc(i)=FINEONLY_X;
	}

	cout << "Ready to build coarse_dog\n";
	/*Dog coarse_dog(coarse_V, coarse_F, coarse_es, coarse_V, coarse_F_ren,
		coarse_submeshVSize, coarse_submeshFSize, submesh_adjacency);
	getInterpolatedCurveCoords(fine_dog,coarse_dog,0);
	return coarse_dog;
	*/
	return Dog(coarse_V, coarse_F, coarse_es, coarse_V, coarse_F_ren,
		coarse_submeshVSize, coarse_submeshFSize, submesh_adjacency);
}

void FineCoarseConversion::print() const {
	/*
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
			cout <<" is "<<ftc[i];
			if(ftc(i) == -2) cout<< " (FINEONLY_I), it lies on coarse edge "<<ftc_edge(i);
			cout<<"\n";
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
	*/
	for(int i=0; i<ftc_curve.size(); ++i){
		for(int j=0;j<ftc_curve[i].size(); ++j){
			cout << "ftc curve ["<<i<<"]["<<j<<"] = "<<ftc_curve[i][j]<<"\n";
		}
	}
}

Eigen::MatrixXd FineCoarseConversion::getCoarseCurveCoords(const Dog& coarse_dog, int curve_idx) const {
	const DogEdgeStitching& es = coarse_dog.getEdgeStitching();
	const Eigen::MatrixXd& V = coarse_dog.getV();
	//Eigen::MatrixXd coarse_coords(entire_coarse_curve_i[curve_idx].size(), 3);
	Eigen::MatrixXd coarse_coords(es.stitched_curves[curve_idx].size(), 3);
	coarse_coords.setZero();
	for(int i=0; i<coarse_coords.rows(); ++i){
		coarse_coords.row(i) = es.stitched_curves[curve_idx][i].getPositionInMesh(V);
		/*
		int row = entire_coarse_curve_i[curve_idx][i];
			for(int j=0; j<4; ++j){
			double weight = entire_coarse_curve_w[curve_idx](row,j);
			if(weight>0)
				coarse_coords.row(i) += weight * V.row(entire_coarse_curve_v[curve_idx](row,j));
		}
		*/
	}
	return coarse_coords;
}

Eigen::MatrixXd FineCoarseConversion::getInterpolatedCurveCoords(const Dog& fine_dog, const Dog& coarse_dog, int curve_idx) const{
	Eigen::MatrixXd coarse_coords = getCoarseCurveCoords(coarse_dog, curve_idx);
	Curve curve(coarse_coords);
	Eigen::RowVector3d T;
	Eigen::Matrix3d F;
	curve.getTranslationAndFrameFromCoords(coarse_coords, T, F);
	Eigen::MatrixXd fine_coords = curve.getInterpolatedCoords(T, F, ctf_curve_offsets[curve_idx]);
	{
		const DogEdgeStitching& es=fine_dog.getEdgeStitching();
		for(int i=0; i<es.stitched_curves[0].size();++i){
			Eigen::RowVector3d p1 = es.stitched_curves[0][i].getPositionInMesh(fine_dog.getV());
			Eigen::RowVector3d p2=fine_coords.row(i);
			double diff = (p1-p2).norm();
			cout << "finest curve 0 "<<i<<": is\t"<<p1<<"\tinterpolated as\t"<<p2;
			if(diff<0.01) cout << "\tgood with diff " << diff <<endl;
			else cout << "\t bad with diff "<<diff<<endl;
		}
	}
	return fine_coords;
}

vector<int> FineCoarseConversion::getFineNeighboursOfCoarseOnly(int coarse) const {
	if(ctf(coarse) >= 0) cout << "Error: Not a coarse-only vertex\n";
	vector<int> res(0);
	//Resulting vector has 0-4 fine LINK vertices, each followed by a FINEONLY_I vertex
	int n = coarseonly_adjacent_links[coarse].size();
	for(int i=0; i<n; ++i){
		res.push_back(coarseonly_adjacent_links[coarse][i]);
		res.push_back(coarseonly_adjacent_fineonly[coarse][i]);
	}
	return res;
}
