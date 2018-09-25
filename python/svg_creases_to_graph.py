import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *

def svg_creases_to_graph(svg_file):
	G = nx.Graph()
	paths, attributes = svg2paths(svg_file)
	print 'Number of paths = ', len(paths)

	viewbox = get_svg_view_box(svg_file)
	print 'viewBox = ', viewbox

	style_classes = get_style_classes(svg_file)
	for i in range(len(paths)):
		path, attrib = paths[i], attributes[i]
		try:
			color = get_curve_color(style_classes,attrib)
			print 'The color is ', color
		except:
			print 'out?'
			continue
		points = get_path_points_as_matrix(path, 5000, float(viewbox[0]), float(viewbox[1]))
		#print 'points = ', points


def test_svg_creases_to_graph():
	svg_creases_to_graph('../crease_patterns/empty.svg')

if __name__ == "__main__":
	test_svg_creases_to_graph()