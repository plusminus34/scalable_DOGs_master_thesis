import networkx as nx
#import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import polygonize,polygonize_full,linemerge,transform
from crease_pattern import *
import svgpathtools.path
import sys

def svg_creases_to_polygonal_data(svg_file):
	G = nx.Graph()
	paths, attributes = svg2paths(svg_file)
	print ('Number of paths = ', len(paths))

	viewbox = get_svg_view_box(svg_file)
	print ('viewBox = ', viewbox)

	style_classes = get_style_classes(svg_file)
	path_lines = []
	for i in range(len(paths)):
		path, attrib = paths[i], attributes[i]
		if is_border(attrib, style_classes):
			print ('Reading border polygon')
			border_poly = handle_border(path)
		else:
			try:
				print ('Reading fold')
				raw_input('am i here')
				vertices = handle_fold(path,5)
				path_lines.append(LineString(vertices))
			except:
				print ('Error handling a fold')

	# Make sure that the border scale is around width = 1, flip the y coordinates (since its opposite in images) translate it such that the lowest coordinate is 0
	border_poly, path_lines = translate_and_normalize_polygons(border_poly, path_lines)
		
	return border_poly,path_lines,viewbox

def handle_border(path):
	if is_polyline(path):
		lines = [((p.start.real, p.start.imag), (p.end.real, p.end.imag)) for p in path]
		"""
		lines = []
		for l in tmp_lines:
			if np.linalg.norm(np.array(l[0]-np.array(l[1]))) > 1e-7:
				print ('adding line = ', l)
				lines.append(l)
			else:
				print ('ignoring l = ', l)
		"""
		#print ('lines = ', lines)
		polys = list(polygonize(lines))
		for p in polys:
			print ('p = ', p)
		assert len(polys) == 1,  'There should be 1 border polygon'
		return polys[0]
	else:
		# bezier curve (not supported for a boundary yet)
		raise NotImplementedError()


# For now assume its a polyline if the first part of the path is (not supporting a path that is a mix of polylines and bezier curves atm)
def is_polyline(path):
	return isinstance(path[0],svgpathtools.path.Line)

def get_border_poly(border_poly):
	convex_hull = MultiPoint(border_poly).convex_hull

def handle_fold(path,sampling,bounds):
	if is_polyline(path):
		print ('Polyline!')
		points = sample_polylines(path,bounds)
		#print ('line points = ', points)
	else:
		print ('bezier curve!')
		print ('sampling = ', sampling)
		points = sample_bezier_path_sampling(path, sampling)
		#print ('bezier points = ', points)
	return points

def translate_and_normalize_polygons(border_poly, path_lines):
	minx, miny, maxx, maxy = border_poly.bounds
	scale = 1./(maxx-minx)
	t = [-scale*minx,scale*maxy]
	#print ('scaling by ', scale)
	
	translate_and_scale = lambda x,y: [x*scale + t[0],-y*scale + t[1]]
	border_poly = transform(translate_and_scale, border_poly) # also invert y coordinates
	new_lines = []
	for p in path_lines:
		new_lines.append(transform(translate_and_scale, p))
	print ('new border_poly.bounds = ', border_poly.bounds)
	return border_poly, new_lines

def test_svg_creases_to_graph(svg_file):
	print ('Testing with file ', svg_file)
	border_poly,polylines = svg_creases_to_polygonal_data(svg_file)
	test_crease_pattern(border_poly, polylines)

if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_svg_creases_to_graph(sys.argv[1])
	else:
		test_svg_creases_to_graph("../crease_patterns/1_curve.svg")
