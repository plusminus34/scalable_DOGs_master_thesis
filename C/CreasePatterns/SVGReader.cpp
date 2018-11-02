#include "SVGReader.h"

#include <igl/readOBJ.h>

#include <sys/stat.h>
#include <stdlib.h>

using namespace std;

void read_svg_crease_pattern(const std::string& path, CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines) {
	// call a python script :P (TODO use a c++ svg reader, will also be way faster)
	const string tmp_folder = "tmp_poly//"; 
	mkdir(tmp_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	const string cmd = string("python ../../python/svg_to_polylines.py ") + path + " " + tmp_folder;
	int res = system(cmd.c_str());
	cout << "cmd returned " << res << endl;

	// Used to read polylines
	Eigen::MatrixXd P; Eigen::MatrixXi Fdummy; 

	// get number of polylines
	int num_poly = 1; 
	while (!ifstream(tmp_folder+to_string(num_poly)).good()) {num_poly++;}
	cout << "num_poly = " << num_poly << endl;
	// Read polylines
	for (int i = 1; i < num_poly+1; i++) {
		igl::readOBJ(tmp_folder+string("poly-")+to_string(i)+string(".obj"),P,Fdummy);
		//cout << "P = " << P << endl;
		auto poly = eigen_to_polyline_2(P);
		//cout << "poly = "<< poly << endl;
		polylines.push_back(poly);
	}

	// read bounding box
	igl::readOBJ(tmp_folder+string("bbox.obj"),P,Fdummy);
	bbox = CGAL::Bbox_2(P(0,0), P(1,0), P(2,0), P(3,0));
	cout << "bbox = " << bbox << endl;

	system(std::string(std::string("rm -r ")+tmp_folder).c_str());
}

Polyline_2 eigen_to_polyline_2(const Eigen::MatrixXd& p) {
	std::vector<Point_2> points(p.rows());
	for (int i = 0; i < p.rows(); i++) {points.push_back(Point_2(p(i,0),p(i,1)));}
	Geom_traits_2 traits; 
	Geom_traits_2::Construct_curve_2 polyline_construct =traits.construct_curve_2_object();
	return polyline_construct(points.begin(),points.end());
}