#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"
#include "../QuadMesh/Quad.h"

class EdgePointConstraints : public Constraints {
public:
	EdgePointConstraints(const QuadTopology& quadTop, std::vector<EdgePoint> edgePoints , const Eigen::MatrixXd& edgePointCoords) :
			vnum(quadTop.v_n), edgePoints(edgePoints) {
		mat2_to_vec(edgePointCoords, bc); // flatten to a vector
		const_n = bc.rows();
		enabled.resize(bc.rows());for(int i=0;i<bc.rows();++i)enabled[i]=true;
		// 2 points per edge
		IJV.resize(2*const_n);
	};
	EdgePointConstraints() {const_n = 0; edgePoints.resize(0); bc.resize(0);enabled.resize(0);} // empty set of constraints

	virtual EdgePointConstraints* clone() const {return new EdgePointConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		int edge_points_n = edgePoints.size(); Eigen::VectorXd coords(3*edge_points_n);
    for (int i = 0; i < edge_points_n; i++) {
      auto vec = edgePoints[i].getPositionInMesh(x,vnum);
      coords(i) = vec(0); coords(edge_points_n+i) = vec(1); coords(2*edge_points_n+i) = vec(2);
    }
		Eigen::VectorXd res = coords-bc;
		for(int i=0;i<enabled.size();++i){
			if(enabled[i])continue;
			res(i)=0;
			res(edge_points_n+i)=0;
			res(2*edge_points_n+i)=0;
		}
    return res;
	}

	void print(const Eigen::VectorXd& x) const {
		int edge_points_n = edgePoints.size(); Eigen::VectorXd coords(3*edge_points_n);
    	for (int i = 0; i < edge_points_n; i++) {
				if(!enabled[i])continue;
        auto vec = edgePoints[i].getPositionInMesh(x,vnum);
        coords(i) = vec(0); coords(edge_points_n+i) = vec(1); coords(2*edge_points_n+i) = vec(2);
				Eigen::RowVector3d target; target<<bc(i),bc(i+edge_points_n),bc(i+2*edge_points_n);
				Eigen::RowVector3d itis; itis<<coords(i),coords(i+edge_points_n),coords(i+2*edge_points_n);
				std::cout << "EdgePointConstraints "<<i<<": should be "<<target<<"\tand is "<<itis<<"\tdiff "<<(target-itis).norm()<<"\n";
    	}
	}

	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_row = 0; int edge_points_n = edgePoints.size(); int ijv_cnt = 0;
		for (int b_i = 0; b_i < edge_points_n; b_i++ ) {
			int v1 = edgePoints[b_i].edge.v1, v2 = edgePoints[b_i].edge.v2;
			double t = edgePoints[b_i].t;
			if(enabled[b_i]){
				// 1 constraint for every coordinate
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_row, v1, t);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_row, v2, 1-t);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(edge_points_n+const_row, vnum+v1, t);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(edge_points_n+const_row, vnum+v2, 1-t);

				IJV[ijv_cnt++] = Eigen::Triplet<double>(2*edge_points_n+const_row, 2*vnum+v1, t);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(2*edge_points_n+const_row, 2*vnum+v2, 1-t);
			}else{
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_row, v1, 0);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(const_row, v2, 0);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(edge_points_n+const_row, vnum+v1, 0);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(edge_points_n+const_row, vnum+v2, 1-t);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(2*edge_points_n+const_row, 2*vnum+v1, 0);
				IJV[ijv_cnt++] = Eigen::Triplet<double>(2*edge_points_n+const_row, 2*vnum+v2, 0);
			}
			const_row++;
		}
	}

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

	void update_coords(Eigen::MatrixXd edgePointCoords) {mat2_to_vec(edgePointCoords, bc);}
	std::vector<EdgePoint> getEdgePoints() const {return edgePoints;}
	Eigen::VectorXd getEdgePointConstraints() const {return bc;}

	void set_enabled(int i, bool val){enabled[i] = val;}

private:
	int vnum;
	std::vector<EdgePoint> edgePoints; Eigen::VectorXd bc;
	std::vector<bool> enabled;// silly little thing to disable curve endpoints for coarse solver
};
