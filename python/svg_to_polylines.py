from svg_utils import *
import sys,os
sys.path.insert(0, os.getcwd() + "/../../libigl/python")
sys.path.insert(0, os.getcwd() + "/../../libigl/external/nanogui/build/python")
import pyigl as igl
from svg_to_polygons import *
from iglhelpers import *


def save_polyline(poly, out_file):
	print 'saving polyline to file ', out_file
	#print 'poly = ', poly
	V = p2e(add_numpy_zero_z_coord_column(poly))
	V.setCol(1,-1*V.col(1))
	#print 'poly = ', V
	igl.writeOBJ(out_file,V,igl.eigen.MatrixXi())

def save_bounding_box(bbox, out_file):
	bbox_mat = np.zeros([4,3])
	bbox_mat[:,0] = bbox
	V = p2e(bbox_mat)
	igl.writeOBJ(out_file,V,igl.eigen.MatrixXi())

def svg_creases_to_bounding_box_and_polylines(svg_file):
	paths, attributes = svg2paths(svg_file)
	print 'Number of paths = ', len(paths)

	viewbox = get_svg_view_box(svg_file)
	print 'viewBox = ', viewbox
	bounds = viewbox[0],viewbox[1],viewbox[0]+viewbox[2],viewbox[1]+viewbox[3]
	print 'bounds = ', bounds

	border_poly = []
	style_classes = get_style_classes(svg_file)
	path_lines = []
	for i in range(len(paths)):
		path, attrib = paths[i], attributes[i]
		print 'path = ', path
		if is_border(attrib, style_classes):
			print 'Reading border polygon'
			border_poly = handle_border(path)
		else:
			try:
				print 'Reading fold'
				vertices = handle_fold(path,500)
				#print 'vertices = ', vertices
				path_lines.append(vertices)
			except:
				print 'Error handling a fold'

	return border_poly,path_lines,bounds

def add_numpy_zero_z_coord_column(V):
	return np.ascontiguousarray(np.vstack((V[:,0],V[:,1],np.zeros(V.shape[0]))).T)

def add_numpy_zero_y_coord_column(V):
	return np.ascontiguousarray(np.vstack((V[:,0],np.zeros(V.shape[0]))).T)

if __name__ == "__main__":
	if len(sys.argv) == 3:
		svg_file, out_folder = sys.argv[1],sys.argv[2]
		border_poly,polylines,boundingBox = svg_creases_to_bounding_box_and_polylines(svg_file)

		print 'boundingBox = ', boundingBox
		save_bounding_box(boundingBox, out_folder+"/bbox.obj")
		#save_polyline(out_folder+"//"+"border_poly.obj")
		cnt = 0
		for poly in polylines:
			save_polyline(poly,out_folder+"/"+"poly-"+str(cnt)+".obj")
	else:
		print 'Usage: svg_to_polylines.py out_folder'