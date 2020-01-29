#include "FineCoarseConversion.h"

#include "../GeometryPrimitives/Curve.h"

FineCoarseConversion::FineCoarseConversion(const Dog& fine_dog, const Dog& coarse_dog){
	const DogEdgeStitching& fine_es = fine_dog.getEdgeStitching();
	const DogEdgeStitching& coarse_es = coarse_dog.getEdgeStitching();
	const Eigen::MatrixXd& fine_V = fine_dog.getV();
	const Eigen::MatrixXi& fine_F = fine_dog.getF();
	const Eigen::MatrixXd& coarse_V = coarse_dog.getV();
	const Eigen::MatrixXi& coarse_F = coarse_dog.getF();
	const QuadTopology& fine_qt = fine_dog.getQuadTopology();
	const QuadTopology& coarse_qt = coarse_dog.getQuadTopology();
	int num_patches = fine_dog.get_submesh_n();

	ftc.resize(fine_V.rows(), -1);
	ctf.resize(coarse_V.rows(), -1);
	ftc_edge.resize(fine_V.rows(), -1);

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
	ftc[fine_origin] = coarse_origin;
	ctf[coarse_origin] = fine_origin;

	// Need to know about vertices near the patch boundary in the flooding
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
	vector<bool> patch_reached(num_patches, false);
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
				ftc[tocheck[i]] = FINEONLY_I;
				//cout << "flood: "<<tocheck[i]<<" is fine-only and has a coarse neighbour\n";
			}
		} else {
			//Even iterations require more thought
			for(int i=0; i<tocheck.size(); ++i){
				ftc[tocheck[i]] = FINEONLY_X;
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
						ctf[coarse_u] = fine_u;
						ftc[fine_u] = coarse_u;
						//cout << "flood: linking coarse "<<coarse_u<<" to fine "<<fine_u <<"\n";
						//Try to cross the border between submeshes
						for(int i=0; i<patch_adj[fine_patch].size(); ++i){
							int other_patch = patch_adj[fine_patch][i];
							if( !patch_reached[other_patch] && fine_vertex_duplicates(fine_u, other_patch) != -1){
								// There are a duplicate of vertices fine_u,coarse_u on other_patch
								patch_reached[other_patch] = true;
								int other_fine_u = fine_vertex_duplicates(fine_u, other_patch);
								int other_coarse_u = coarse_vertex_duplicates(coarse_u, other_patch);
								ftc[other_fine_u] = other_coarse_u;
								ctf[other_coarse_u] = other_fine_u;
								// Check neighbours on the other patch in the next iteration
								for(int j=0; j<fine_vv[other_fine_u].size(); ++j){
									int willcheck = fine_vv[other_fine_u][j];
									ftc[willcheck] = WILLBECHECKED;
									tocheck_next.push_back(willcheck);
								}
								for(int j=0; j<coarse_vv[other_coarse_u].size(); ++j){
									int willcheck = coarse_vv[other_coarse_u][j];
									ctf[willcheck] = WILLBECHECKED;
//									coarse_tocheck_next.push_back(willcheck);
								}
							}
						}

					}
				}
				/*
				for(int j=0; j<coarse_vv[coarse_tocheck[coarse_i]].size(); ++j){
					int coarse_u = coarse_vv[coarse_tocheck[coarse_i]][j];
					if(ctf[coarse_u] == UNCHECKED){
						ctf[coarse_u] = WILLBECHECKED;
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
				if(ftc[v] == UNCHECKED){
					ftc[v] = WILLBECHECKED;
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

	// get bounding boxes of patches
	vector<double> min_x(num_patches), max_x(num_patches), min_y(num_patches), max_y(num_patches);
	double dx = 0.0, dy = 0.0;
	for(int i=0; i<num_patches; ++i){
		int mini,maxi;
		coarse_dog.get_submesh_min_max_i(i, mini,maxi);
		min_x[i] = max_x[i] = coarse_V(mini, 0);
		min_y[i] = max_y[i] = coarse_V(mini, 1);
		for(int j = mini + 1; j <= maxi; ++j){
			if(coarse_V(j,0) < min_x[i]) min_x[i] = coarse_V(j,0);
			else if(coarse_V(j,0) > max_x[i]) max_x[i] = coarse_V(j,0);
			if(coarse_V(j,1) < min_y[i]) min_y[i] = coarse_V(j,1);
			else if(coarse_V(j,1) > max_y[i]) max_y[i] = coarse_V(j,1);
		}
	}
	// Assumption: Distances between adjacent vertices are constant across the Dog
	for(int i=0; i<fine_vv[0].size(); ++i){
		dx = max(dx, abs(fine_V(0,0) - fine_V(fine_vv[0][i],0)) );
		dy = max(dy, abs(fine_V(0,1) - fine_V(fine_vv[0][i],1)) );
	}
	// construct grids
	vector<Eigen::MatrixXi> fine_grid(num_patches);
	vector<Eigen::MatrixXi> coarse_grid(num_patches);
	for(int i=0; i<num_patches; ++i){
		//create grid
		int fine_x_dim = round( (max_x[i] - min_x[i]) / dx ) + 1;
		int fine_y_dim = round( (max_y[i] - min_y[i]) / dy ) + 1;
		fine_grid[i].setConstant(fine_x_dim, fine_y_dim, -1);
		int coarse_x_dim = ceil(0.5 * fine_x_dim);
		int coarse_y_dim = ceil(0.5 * fine_y_dim);
		coarse_grid[i].setConstant( coarse_x_dim, coarse_y_dim, -1);

		//fill vertices into grid
		int mini,maxi;
		fine_dog.get_submesh_min_max_i(i, mini,maxi);
		for(int j = mini; j <= maxi; ++j){
			int x_idx = round( (fine_V(j,0)-min_x[i])/dx );
			int y_idx = round( (fine_V(j,1)-min_y[i])/dy );
			fine_grid[i](x_idx, y_idx) = j;
		}
		coarse_dog.get_submesh_min_max_i(i, mini,maxi);
		for(int j = mini; j <= maxi; ++j){
			int x_idx = round( (coarse_V(j,0)-min_x[i]) / (2*dx) );
			int y_idx = round( (coarse_V(j,1)-min_y[i]) / (2*dy) );
			coarse_grid[i](x_idx, y_idx) = j;
		}
		//mark places where curve submeshes need additional vertices
		//(Of course, curve submeshes got removed. But the labels are still useful.)
		for(int j = mini; j <= maxi; ++j){
			if(ctf[j] > -1 ) continue;// coarse-only vertices are relevant
			int x_idx = round( (coarse_V(j,0)-min_x[i]) / (2*dx) );
			int y_idx = round( (coarse_V(j,1)-min_y[i]) / (2*dy) );
			fine_grid[i](2*x_idx, 2*y_idx) = CURVEONLY_LINK;
			//mark fine vertices between adjacent coarse vertices
			if(x_idx > 0 && coarse_grid[i](x_idx-1, y_idx) > -1 && fine_grid[i](2*x_idx-1, 2*y_idx) == -1)
				fine_grid[i](2*x_idx-1, 2*y_idx) = CURVEONLY_I;
			if(y_idx > 0 && coarse_grid[i](x_idx, y_idx-1) > -1 && fine_grid[i](2*x_idx, 2*y_idx-1) == -1)
				fine_grid[i](2*x_idx, 2*y_idx-1) = CURVEONLY_I;
			if(x_idx < coarse_x_dim-1 && coarse_grid[i](x_idx+1, y_idx) > -1 && fine_grid[i](2*x_idx+1, 2*y_idx) == -1)
				fine_grid[i](2*x_idx+1, 2*y_idx) = CURVEONLY_I;
			if(y_idx < coarse_y_dim-1 && coarse_grid[i](x_idx, y_idx+1) > -1 && fine_grid[i](2*x_idx, 2*y_idx+1) == -1)
				fine_grid[i](2*x_idx, 2*y_idx+1) = CURVEONLY_I;
		}
	}
	// use grids to map fine I vertices to coarse edges
	for(int i=0; i<coarse_qt.E.rows(); ++i){
		double x0 = coarse_V(coarse_qt.E(i,0), 0);
		double y0 = coarse_V(coarse_qt.E(i,0), 1);
		double x1 = coarse_V(coarse_qt.E(i,1), 0);
		double y1 = coarse_V(coarse_qt.E(i,1), 1);
		int patch = coarse_dog.v_to_submesh_idx(coarse_qt.E(i,0));
		int x0_idx = round( (x0-min_x[patch]) / (2*dx) );
		int y0_idx = round( (y0-min_y[patch]) / (2*dy) );
		int x1_idx = round( (x1-min_x[patch]) / (2*dx) );
		int y1_idx = round( (y1-min_y[patch]) / (2*dy) );
		// idx in fine_grid is just the sum of the coarse idx!
		int fine = fine_grid[patch](x0_idx + x1_idx, y0_idx + y1_idx);
		if(fine > -1 && fine < ftc_edge.size()) ftc_edge[fine] = i;
	}

	// prepare data for fine-to-coarse update
	ftc_update_vertices.resize(ctf.size());
	ftc_update_weights.resize(ctf.size());
	for(int i=0; i<ctf.size(); ++i){
		if(ctf[i] < 0){
			// coarse only vertices have to be approximated
			int patch = coarse_dog.v_to_submesh_idx(i);
			int x_idx = round( (coarse_V(i,0) - min_x[patch]) / dx );//in fine grid
			int y_idx = round( (coarse_V(i,1) - min_y[patch]) / dy );
			cout << "coarse vertex "<<i<<" is at finegrid "<<x_idx<<" "<<y_idx<<endl;
			cout<<"    coarse row "<<i<<"   "<<coarse_V.row(i)<<endl;
			int num_fine_quads = 0;
			// check adjacent coarse quads;
			for(int j=0; j<coarse_qt.VF[i].size(); ++j){
				int cq = coarse_qt.VF[i][j];// row cq of coarse_F
				//cout << "  coarse quad "<<cq <<" is adjacent\n";
				// get the grid coordinates of the associated fineonly-X vertex
				double x_mid = 0.5*(coarse_V(coarse_F(cq,0), 0) + coarse_V(coarse_F(cq,2), 0));
				double y_mid = 0.5*(coarse_V(coarse_F(cq,0), 1) + coarse_V(coarse_F(cq,2), 1));
				int xp_x_idx = round( (x_mid - min_x[patch]) / dx );
				int xp_y_idx = round( (y_mid - min_y[patch]) / dy );
				int xp = fine_grid[patch](xp_x_idx, xp_y_idx);
				int k;
				int close_links = 0;
				for(int l=0;l<4;++l){
					if(coarse_F(cq,l) == i) {k = l;}
					if(ctf[coarse_F(cq,l)] > -1) ++close_links;
				}
				if(close_links==3){
					// do it with the close link points
					++num_fine_quads;
					int opposite = ctf[coarse_F(cq, (k+2)%4)];
					int u1 = ctf[coarse_F(cq, (k+1)%4)];
					int u2 = ctf[coarse_F(cq, (k+3)%4)];
					ftc_update_vertices[i].push_back(opposite);
					ftc_update_weights[i].push_back(-1.0);
					ftc_update_vertices[i].push_back(u1);
					ftc_update_weights[i].push_back(1.0);
					ftc_update_vertices[i].push_back(u2);
					ftc_update_weights[i].push_back(1.0);
				} else {
					for(int l=1; l<4; ++l){
						int ll = (k+l) % 4;
						if(ctf[coarse_F(cq, ll)] < 0) continue;
						// LINK-vertex implies a fine quad within the coarse quad
						//   (except that's not quite right, see noncorner case)
						++num_fine_quads;
						int link_vertex = ctf[coarse_F(cq, ll)];
						int link_x_idx = round( (fine_V(link_vertex, 0) - min_x[patch]) / dx );
						int link_y_idx = round( (fine_V(link_vertex, 1) - min_y[patch]) / dy );
						if(l == 2){
							int other_a = fine_grid[patch](xp_x_idx, link_y_idx);
							int other_b = fine_grid[patch](link_x_idx, xp_y_idx);
							if(other_a > -1 && other_b > -1){
								//it's a corner
								// variant "minimal vertices"
								ftc_update_vertices[i].push_back(link_vertex);
								ftc_update_weights[i].push_back(-1.0);
								ftc_update_vertices[i].push_back(xp);
								ftc_update_weights[i].push_back(2.0);
								//TODO variant 2
							} else {
								--num_fine_quads;
							}
						} else {
							int midpoint = fine_grid[patch]((x_idx+link_x_idx)/2, (y_idx+link_y_idx)/2);
							if(midpoint>-1){
								//not a corner
								// variant "minimal vertices"
								ftc_update_vertices[i].push_back(link_vertex);
								ftc_update_weights[i].push_back(-1.0);
								ftc_update_vertices[i].push_back(midpoint);
								ftc_update_weights[i].push_back(2.0);
							} else {
								--num_fine_quads;
							}
						}
					}
				}
			}
			for(int j=0; j<ftc_update_weights[i].size(); ++j){
				// equal weights for all close fine quads
				if(num_fine_quads != 0) ftc_update_weights[i][j] /= num_fine_quads;
			}
		} else {
			// link vertices are easy
			cout << "coarse vertex "<<i<<" is a link point\n";
			//cout <<"pushback link "<<ctf[i]<<endl;
			ftc_update_vertices[i].push_back(ctf[i]);
			ftc_update_weights[i].push_back(1);
		}
	}
	/*
	for(int i=0; i<num_patches;++i){
		cout<<"fine patch "<<i<<":\n"<<fine_grid[i]<<endl;
		cout<<"coarse patch "<<i<<":\n"<<coarse_grid[i]<<endl;
	}
	*/
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
			if(ftc[i] == -2) cout<< " (FINEONLY_I), it lies on coarse edge "<<ftc_edge[i];
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
	return fine_coords;
}

Eigen::MatrixXd FineCoarseConversion::coarsen(const Eigen::MatrixXd& fine_V) const{
	Eigen::MatrixXd coarse_V;
	coarse_V.setZero(ctf.size(), 3);
	for(int i=0; i<coarse_V.rows(); ++i){
		for(int j=0; j<ftc_update_weights[i].size(); ++j){
			coarse_V.row(i) += ftc_update_weights[i][j] * fine_V.row(ftc_update_vertices[i][j]);
		}
	}
	return coarse_V *0.5;//coarse scale
}
