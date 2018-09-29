import matplotlib.pyplot as plt
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import polygonize,polygonize_full,linemerge,transform
from crease_pattern import *
from svg_to_polygons import *
from planar_dog import *
import sys

# need also polyline since part of the polygons is actually the face border, making it all a bit more complicated, I think
def dog_from_border_and_polylines(border_poly, polylines):
	pass

def compute_grid_border_polygon(border_poly, grid):
	pass

def test_dog_from_face_polygons(svg_file):
	print 'Testing with file ', svg_file
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

	border_poly,polylines = svg_creases_to_polygonal_data(svg_file)
	face_polygons, polylines = crease_pattern(border_poly, polylines)

	res_x,res_y = 25,25
	grid = grid_from_boundary(border_poly, res_x,res_y)
	plot_face_polygons(face_polygons, polylines, ax1, 'Faces with grid')
	plot_grid(grid, ax1)

	dog = dog_from_border_and_polylines(border_poly, polylines)

	# show all
	plt.show()


if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_from_face_polygons(sys.argv[1])
	else:
		test_dog_from_face_polygons("../crease_patterns/1_curve.svg")
