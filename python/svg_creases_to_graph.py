import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import polygonize,polygonize_full,linemerge
from crease_pattern import *
import svgpathtools.path
from sys import exit

def svg_creases_to_polygonal_data(svg_file):
	G = nx.Graph()
	paths, attributes = svg2paths(svg_file)
	print 'Number of paths = ', len(paths)

	viewbox = get_svg_view_box(svg_file)
	print 'viewBox = ', viewbox

	style_classes = get_style_classes(svg_file)
	path_lines = []
	for i in range(len(paths)):
		path, attrib = paths[i], attributes[i]
		if is_border(attrib, style_classes):
			border_poly = handle_border(path)
		else:
			try:
				vertices = handle_fold(path,100)
				path_lines.append(LineString(vertices))
			except:
				print 'Error handling a fold'
		
	return border_poly,path_lines

def handle_border(path):
	if is_polyline(path):
		lines = [((p.start.real, p.start.imag), (p.end.real, p.end.imag)) for p in path]
		#print 'lines = ', lines
		polys = list(polygonize(lines))
		assert len(polys) == 1,  'There should be 1 border polygon'
		return polys[0]
	else:
		# bezier curve
		pass
def is_border(attrib,style_classes):
	color = get_curve_color(style_classes,attrib)
	print 'The color is ', color
	is_border = (color == (0,0,0))
	return is_border

# For now assume its a polyline if the first part of the path is (not supporting a path that is a mix of polylines and bezier curves atm)
def is_polyline(path):
	return isinstance(path[0],svgpathtools.path.Line)

def get_border_poly(border_poly):
	#print 'border_poly =', border_poly
	#bla = shapely.convex_hull(border_poly)

	#bla = LinearRing([(0, 0), (1, 1), (1, 0)])
	convex_hull = MultiPoint(border_poly).convex_hull

def handle_fold(path,sampling = 500):
	if is_polyline(path):
		print 'this is a line!'
		points = sample_polylines(path)
		print 'points = ', points
		exit(1)
	else:
		print 'bezier curve!'
		points = sample_bezier_path_sampling(path, sampling)
		print 'points.shape = ', points.shape
		#fdsfd
	return points

def test_svg_creases_to_graph():
	
	border_poly,polylines = svg_creases_to_polygonal_data('../crease_patterns/empty.svg')
	#test_plot_polygon_and_lines(1,border_poly,[polylines])	
	# hardcoded..
	#viewBox =  [0.0, -1000.0, 1500.0, 1000.0]
	#border_poly = Polygon([(132,868),(1370,868),(1370,153),(132,153)])
	test_crease_pattern(border_poly, polylines)


if __name__ == "__main__":
	test_svg_creases_to_graph()