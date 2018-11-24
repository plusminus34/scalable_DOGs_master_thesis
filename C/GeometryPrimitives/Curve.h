#pragma once

#include <Eigen/Dense>
#include "../QuadMesh/Quad.h"

class SurfaceCurve {
public:
	Eigen::MatrixXd get_curve_coords(const Eigen::MatrixXd& V) const {
		Eigen::MatrixXd coords(edgePoints.size(),3);
		for (int i = 0; i < edgePoints.size(); i++) {
			coords.row(i) = edgePoints[i].getPositionInMesh(V);
		}
		return coords;
	}
	std::vector<EdgePoint> edgePoints;
};

// Describe a curve up to rigid motion with discrete edge lengths, curvature and torsion
class Curve {
public:
	Curve(const std::vector<double>& len, const std::vector<double>& k, const std::vector<double>& t);
	Curve(const Eigen::MatrixXd& coords);
	Curve(const SurfaceCurve& surfaceCurve, const Eigen::MatrixXd& V);
	// Interpolate curves
	Curve(const std::vector<double>& len1, const std::vector<double>& k1, const std::vector<double>& t1,
		const std::vector<double>& len2, const std::vector<double>& k2, const std::vector<double>& t2, double time);

	Curve(const Curve& curve1, const Curve& curve2, double time);

	// Return a curve coordinates given rigid motion data (translation and rotation frame)
	// The translation and rotation are for the second point on the curve (need a point with an osculation plane)
	Eigen::MatrixXd getCoords(const Eigen::RowVector3d& T, const Eigen::Matrix3d& F);

	static void getTranslationAndFrameFromCoords(const Eigen::MatrixXd& coords, Eigen::RowVector3d& T, Eigen::Matrix3d& F);
	
	// edge lengths, vertex curvature, edge torsion
	std::vector<double> len, k, t;
private:
	double get_angle(Eigen::RowVector3d e1,Eigen::RowVector3d e2);
	double get_angle_from_lengths_and_k(double l1, double l2, double k);


	double euc_to_arc(double arc, double k);
	double arc_to_euc(double s, double k);
};
