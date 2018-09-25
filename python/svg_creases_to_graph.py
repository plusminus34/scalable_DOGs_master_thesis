import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg

def svg_creases_to_graph(svg_file):
	G = nx.Graph()
	paths, attributes = svg2paths(svg_file)

def test_svg_creases_to_graph():
	pass

if __name__ == "__main__":
	test_svg_creases_to_graph()