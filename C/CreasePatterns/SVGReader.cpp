#include "SVGReader.h"

#include <igl/readOBJ.h>

#include <sys/stat.h>
#include <stdlib.h>

using namespace std;

void read_svg_crease_pattern(const std::string& path, CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines,
		std::vector<Polyline_2>& boundary_polylines, std::vector<Eigen::MatrixXd>& polylines_data,
		std::vector<Eigen::MatrixXd>& boundary_polylines_data) {
	// call a python script :P (TODO use a c++ svg reader, will also be way faster)
	const string tmp_folder = "tmp_poly//";
	mkdir(tmp_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	const string cmd = string("python ../../python/svg_to_polylines.py ") + path + " " + tmp_folder;
	int res = system(cmd.c_str());
	cout << "cmd returned " << res << endl;

	// Read polylines
	polylines_data.clear();
	Eigen::MatrixXd P; Eigen::MatrixXi Fdummy;
	int next_poly = 1;
	while (igl::readOBJ(tmp_folder+string("poly-")+to_string(next_poly)+string(".obj"),P,Fdummy)) {
		auto poly = points_to_polylines(P);
		polylines.push_back(poly);
		polylines_data.push_back(P);
		next_poly++;
	}

	boundary_polylines_data.clear();
	int next_border_poly = 1;
	while (igl::readOBJ(tmp_folder+string("borderpoly-")+to_string(next_border_poly)+string(".obj"),P,Fdummy)) {
		auto poly = points_to_polylines(P);
		boundary_polylines.push_back(poly);
		boundary_polylines_data.push_back(P);
		next_border_poly++;
	}

	if (boundary_polylines.size() == 0) {
		// Read bounding box
		igl::readOBJ(tmp_folder+string("bbox.obj"),P,Fdummy);
		bbox = CGAL::Bbox_2(P(0,0), P(1,0), P(2,0), P(3,0));
		cout << "bbox = " << bbox << endl;

		// add a boundary polyline using the boundary box
		double x0 = bbox.xmin(), x1 = bbox.xmax(), y0 = bbox.ymin(), y1 = bbox.ymax();
		Point_2 pt1(x0,y0), pt2(x0,y1), pt3(x1, y1), pt4(x1,y0);
		std::list<Point_2> pts = {pt1,pt2,pt3,pt4,pt1}; // circular list
		Geom_traits_2 traits;
		Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();
		boundary_polylines.push_back(polyline_construct(pts.begin(), pts.end()));
		boundary_polylines_data.resize(1);
		Eigen::MatrixXd P_data(5,3);
		P_data << x0, y0, 0,
							x0, y1, 0,
							x1, y1, 0,
							x1, y0, 0,
							x0, y0, 0;
		boundary_polylines_data[0] = P_data;
	} else {
		// Set boundary box from boundary polylines
		bbox = boundary_polylines[0].bbox();
		for (int i = 1; i < boundary_polylines.size(); i++) bbox += boundary_polylines[i].bbox();
	}
	system(std::string(std::string("rm -r ")+tmp_folder).c_str());
}

//Polyline_2 points_to_polylines_snapped_at_start_end(const Eigen::MatrixXd& p, const CGAL::Bbox_2& bbox) {
Polyline_2 points_to_polylines(const Eigen::MatrixXd& p) {
	std::vector<Point_2> points(p.rows());
	for (int i = 0; i < p.rows(); i++) {points[i] = Point_2(p(i,0),p(i,1));}
	Geom_traits_2 traits;
	Geom_traits_2::Construct_curve_2 polyline_construct =traits.construct_curve_2_object();

	int l_idx = p.rows()-1;


	//points[0] = snap_pt_to_bbox(points[0], bbox);
	//points[l_idx] = snap_pt_to_bbox(points[l_idx], bbox);

	Polyline_2 poly =  polyline_construct(points.begin(),points.end());
	/*
	std::cout << "printing poly with poly points: " << std::endl; for (auto v: points) std::cout << "v = " << v << std::endl;
	std::cout << "printing poly with poly: " << std::endl; std::cout << "v = " << poly.subcurves_begin()->source() << std::endl;
	for (auto it = poly.subcurves_begin(); it != poly.subcurves_end(); it++) std::cout << it->target() << std::endl;
	int wait; cin >> wait;
	*/
	return poly;
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
