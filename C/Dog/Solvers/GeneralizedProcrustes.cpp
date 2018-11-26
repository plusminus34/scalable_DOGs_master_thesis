#include "GeneralizedProcrustes.h"

#include <igl/procrustes.h>

double GeneralizedProcrustes::solve(Dog& dog, int fixed_mesh_i,
		const PositionalConstraints& posConst,
        const StitchingConstraints& stitchingConstraints,
        /*const PositionalConstraints& EdgePointConstraints, TODO */
        Eigen::MatrixXd& Vout) {

	procrustes_on_submesh(dog, fixed_mesh_i, posConst);

	// For now just support 2 meshes, so we assume all the linear constraints are determined
	// In general we need to do a bfs, each time propagating constraints and setting free variables after procrustes
	// Saving which free variables we already have set
	int submesh_i = dog.get_submesh_n();
	for (int i = 0; i < submesh_i; i++) {
		if (i!=fixed_mesh_i) {
			// For the second mesh we have a lot of stitching constraints,
			// Since we fixed the other submesh these are now edgePoint constraints
			// So we have edgePoint constraints, as well as a point constraint
			// We can now do regular procrustes on the substitution of variables

			// This means we need a function to handle edgePoint constraints together with positional constraints
			//  on a submesh (should be almost the same as procrustes_on_submesh, so we can actually just make it bigger)


			// For now assume only 2 meshes (so all stitching constraints can be changed to edge point ones)
			const DogEdgeStitching& eS = dog.getEdgeStitching();
			std::vector<EdgePoint> edgePoints(eS.edge_const_1.size()); Eigen::MatrixXd edgePointCoords(edgePoints.size(),3);
			//EdgePointConstraints(std::vector<EdgePoint> edgePoints , const Eigen::MatrixXd& edgePointCoords)
			for (int i = 0; i < eS.edge_const_1.size(); i++) {
				//if (eS.edge_const_1[i].v1)
				//std::vector<Edge> edge_const_1, edge_const_2;
			}
	//std::vector<double> edge_coordinates;
		}
	}
	return 0;
}

void GeneralizedProcrustes::procrustes_on_submesh(Dog& dog, int submesh_i, 
													const PositionalConstraints& posConst,
													const EdgePointConstraints& edgePointConstraints) {
	Eigen::MatrixXd& dogV = dog.getVMutable();
	// Build set of positional constraint points in the mesh
	Eigen::VectorXi b(posConst.getPositionIndices()); Eigen::VectorXd bc(posConst.getPositionVals());
	Eigen::MatrixXd targetPointPositions; vec_to_mat2(bc,targetPointPositions);
	Eigen::VectorXi bV(b.rows()/3); for (int i = 0; i < bV.rows(); i++) bV(i) = b(i);
	Eigen::MatrixXd constrainedPositions; igl::slice(dog.getV(),bV,1, constrainedPositions);
	
	
	// Build set of positional constraint points from edge point constraints
	Eigen::VectorXd edgePointsBc(edgePointConstraints.getEdgePointConstraints());
	Eigen::MatrixXd targetEdgePoints; vec_to_mat2(edgePointsBc,targetEdgePoints);
	std::vector<EdgePoint> edgePoints = edgePointConstraints.getEdgePoints();
	Eigen::MatrixXd constrainedEdgePoints = EdgePoint::getPositionInMesh(edgePoints, dog.getV());


	Eigen::MatrixXd src(constrainedPositions.rows() + constrainedEdgePoints.rows(),3);
	src << constrainedPositions,constrainedEdgePoints;
	Eigen::MatrixXd target(targetPointPositions.rows()+targetEdgePoints.rows(), 3);
	target << targetPointPositions, targetEdgePoints;


	// rigid motion (and scale which we ignore for now)
	Eigen::MatrixXd R; Eigen::VectorXd t; double scale_dummy;
	igl::procrustes(src,target,false,false,scale_dummy,R,t);

	// go through mesh vertices and set them
	int submesh_min_i, submesh_max_i; dog.get_submesh_min_max_i(submesh_i, submesh_min_i, submesh_max_i);
	auto Rt = R.transpose();
	for (int i = submesh_min_i; i <= submesh_max_i; i++) {
		// Todo check if we need Rt or R here..c
		dogV.row(i) = (dogV.row(i)*R)+t.transpose();
	}
}

void GeneralizedProcrustes::procrustes_on_submesh(Dog& dog, int submesh_i, 
													const PositionalConstraints& posConst) {
	EdgePointConstraints edgePointConst; // empty edge point constraints
	procrustes_on_submesh(dog, submesh_i, posConst, edgePointConst);
}