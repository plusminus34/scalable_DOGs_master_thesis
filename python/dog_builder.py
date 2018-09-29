import matplotlib.pyplot as plt
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import polygonize,polygonize_full,linemerge,transform
from crease_pattern import *
from svg_to_polygons import *
from planar_dog import *
import sys

# need also polyline since part of the polygons is actually the face border, making it all a bit more complicated, I think
def intersected_grid_and_polylines(grid, polylines):
	grid_intersected_polylines = []
	poly_intersections = find_polylines_intersections(polylines)
	grid = split_grid_by_intersections(grid, poly_intersections)
	for p in polylines:
		#print 'p.intersects(grid) = ', p.intersects(GeometryCollection(grid))
		grid_int = p.intersection(GeometryCollection(grid))
		grid_intersected_polylines.append(LineString(grid_int))
		#print 'grid int len = ', len(list(grid_int))
	return grid, grid_intersected_polylines

def find_polylines_intersections(polylines):
	int_points = []
	for pol_line in polylines:
		new_line = pol_line
		for pol_line2 in polylines:
			if pol_line != pol_line2:
				lines_int = pol_line.intersection(pol_line2)
				if isinstance(lines_int,GeometryCollection):
					for p in lines_int:
						int_points = int_points + p.coords[:]
				else:
					int_points = int_points + lines_int.coords[:]
	# get unique vertices
	int_points = unique_rows(int_points)
	#print 'int_points = ', int_points
	return int_points


def build_mesh_from_grid_and_polylines(grid, polylines):
	pass

def compute_grid_intersected_polygons(grid, border_poly, polylines):
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

	#dog = dog_from_border_and_polylines(border_poly, polylines)
	grid,grid_polylines = intersected_grid_and_polylines(grid, polylines)
	plot_face_polygons(face_polygons, polylines, ax2, 'Grid intersections')
	plot_grid(grid, ax2)
	for line in grid_polylines:
		plot_coords(ax2, line)
		plot_line(ax2, line, 1, '#ffffff') # line width = 1


	border_poly,grid_polylines = remove_points_outside_border(border_poly, grid_polylines)
	face_polygons = build_polygons(border_poly, grid_polylines)
	face_polygons_num = len(list(face_polygons))
	plot_face_polygons(face_polygons, grid_polylines, ax3, 'Faces polyline decomposition (' + str(face_polygons_num) + ' faces)')

	# show all
	plt.show()


if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_from_face_polygons(sys.argv[1])
	else:
		test_dog_from_face_polygons("../crease_patterns/1_curve.svg")
