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

	// Used to read the bounding box and polylines
	Eigen::MatrixXd P; Eigen::MatrixXi Fdummy; 

	// read bounding box
	igl::readOBJ(tmp_folder+string("bbox.obj"),P,Fdummy);
	bbox = CGAL::Bbox_2(P(0,0), P(1,0), P(2,0), P(3,0));
	cout << "bbox = " << bbox << endl;

	// get number of polylines
	int num_poly = 1; 
	while (!ifstream(tmp_folder+to_string(num_poly)).good()) {num_poly++;}
	cout << "num_poly = " << num_poly << endl;
	// Read polylines
	for (int i = 1; i < num_poly+1; i++) {
		igl::readOBJ(tmp_folder+string("poly-")+to_string(i)+string(".obj"),P,Fdummy);
		//cout << "P = " << P << endl;
		auto poly = points_to_polylines_snapped_at_start_end(P,bbox);
		//cout << "poly = "<< poly << endl;
		polylines.push_back(poly);
	}

	//system(std::string(std::string("rm -r ")+tmp_folder).c_str());
}

Polyline_2 points_to_polylines_snapped_at_start_end(const Eigen::MatrixXd& p, const CGAL::Bbox_2& bbox) {
	std::vector<Point_2> points(p.rows());
	for (int i = 0; i < p.rows(); i++) {points[i] = Point_2(p(i,0),p(i,1));}
	Geom_traits_2 traits; 
	Geom_traits_2::Construct_curve_2 polyline_construct =traits.construct_curve_2_object();

	int l_idx = p.rows()-1;

	
	points[0] = snap_pt_to_bbox(points[0], bbox);
	points[l_idx] = snap_pt_to_bbox(points[l_idx], bbox);
	
	return polyline_construct(points.begin(),points.end());
}

Point_2 snap_pt_to_bbox(const Point_2& pt, const CGAL::Bbox_2& bbox) {
	Number_type x = pt.x(),y = pt.y();
	//cout << endl << endl << "------pt = " << pt << "-------"<< endl <<endl;
	round_if_close(x,bbox.xmin());
	round_if_close(x,bbox.xmax());
	round_if_close(y,bbox.ymin());
	round_if_close(y,bbox.ymax());

	return Point_2(x,y);
}

void round_if_close(Number_type& pt, const Number_type& floor_pt) {
	Number_type eps(1e-10);
	bool is_close = (CGAL::abs(pt-floor_pt) < eps);
	if ((is_close) && (pt!=floor_pt)) {
		std::cout << "snapping " << pt << " to " << floor_pt << endl;
		pt = floor_pt;
	}
}