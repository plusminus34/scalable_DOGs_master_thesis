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
		int_coords = p.intersection(GeometryCollection(grid))
		grid_int = np.array([pt.coords[0] for pt in int_coords])
		grid_int = sort_grid_int_by_polyline_points(p, grid_int)
		grid_intersected_polylines.append(LineString(grid_int))
		#print 'grid int len = ', len(list(grid_int))
	return grid, grid_intersected_polylines

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

def sort_grid_int_by_polyline_points(polyline, grid_int):
	grid_int_sorted_coords = []
	pt_sample = 5*grid_int.shape[0]
	poly_coords = [polyline.interpolate(t,True).coords[:] for t in np.linspace(0,1,pt_sample
		)]
	#polyline.coords[:]

	
	distances_idx = np.empty(len(poly_coords))
	i = 0
	for pt in poly_coords:
		#print 'pt = ', pt
		distances_idx[i] = closest_node(pt, grid_int)
		i += 1
	#print 'distances_idx = ', distances_idx
	cur_idx = distances_idx[0]
	
	pt = grid_int[int(cur_idx)]
	cnt = 1
	#grid_int.remove(pt)
	grid_int_sorted_coords.append(pt)
	for idx in distances_idx:
		if idx != cur_idx:
			cur_idx = int(idx)
			pt = grid_int[int(cur_idx)]
			#indices = np.where()
			#grid_int.remove(pt)
			grid_int_sorted_coords.append(pt)
			cnt+=1

	assert len(grid_int_sorted_coords) == grid_int.shape[0], 'Error: lengths do not match'
	return grid_int_sorted_coords

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
	plot_grid(grid, ax1,1.5, '#ffffff')

	#dog = dog_from_border_and_polylines(border_poly, polylines)
	grid,grid_polylines = intersected_grid_and_polylines(grid, polylines)
	plot_face_polygons(face_polygons, polylines, ax2, 'Grid intersections')
	plot_grid(grid, ax2, 1.5, '#ffffff')
	for line in grid_polylines:
		plot_coords(ax2, line, '#bbbbff')
		plot_line(ax2, line, 1, '#dddddd') # line width = 1


	border_poly,grid_polylines = remove_points_outside_border(border_poly, grid_polylines)
	face_polygons = build_polygons(border_poly, grid_polylines)
	face_polygons_num = len(list(face_polygons))
	plot_face_polygons(face_polygons, grid_polylines, ax3, 'Faces polyline decomposition (' + str(face_polygons_num) + ' faces)')
	#pol_patch = PolygonPatch(border_poly)
	#ax3.add_patch(pol_patch)
	plot_line(ax3,LineString(border_poly.exterior.coords),1,'#000000')

	# show all
	plt.show()


if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_from_face_polygons(sys.argv[1])
	else:
		test_dog_from_face_polygons("../crease_patterns/1_curve.svg")
