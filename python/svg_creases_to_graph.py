import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *
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
		#print 'path number ', i, ' is ', path, ' dir(path) = ', dir(path)
		try:
			vertices, is_border = handle_path(path,attrib,style_classes,viewbox,100)
			print 'is_border = ', is_border
			if is_border:
				border_poly = Polygon(vertices) # Take the border as a polygon
			else:
				bla = LineString(vertices)
				path_lines.append(bla)
				pass
		except:
			print 'error handling one path'
	return border_poly,path_lines


def get_border_poly(border_poly):
	#print 'border_poly =', border_poly
	#bla = shapely.convex_hull(border_poly)

	#bla = LinearRing([(0, 0), (1, 1), (1, 0)])
	convex_hull = MultiPoint(border_poly).convex_hull

def handle_path(path,attrib,style_classes,viewbox,sampling = 500):
	color = get_curve_color(style_classes,attrib)
	print 'The color is ', color
	if isinstance(path[0],svgpathtools.path.Line):
		print 'this is a line!'
		points = sample_polylines(path)
		print 'points = ', points
		exit(1)
	else:
		print 'bezier curve!'
		points = sample_bezier_path_sampling(path, sampling)
		print 'points.shape = ', points.shape
		#fdsfd
	is_border = (color == (0,0,0))
	return points,is_border

def test_svg_creases_to_graph():
	
	border_poly,polylines = svg_creases_to_polygonal_data('../crease_patterns/empty.svg')
	#test_plot_polygon_and_lines(1,border_poly,[polylines])	
	# hardcoded..
	#viewBox =  [0.0, -1000.0, 1500.0, 1000.0]
	border_poly = Polygon([(132,868),(1370,868),(1370,153),(132,153)])
	test_crease_pattern(border_poly, polylines)


if __name__ == "__main__":
	test_svg_creases_to_graph()