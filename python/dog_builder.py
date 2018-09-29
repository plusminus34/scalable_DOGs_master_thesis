import matplotlib.pyplot as plt
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import polygonize,polygonize_full,linemerge,transform
from crease_pattern import *
from svg_to_polygons import *
import sys

# need also polyline since part of the polygons is actually the face border, making it all a bit more complicated, I think
def dog_from_face_polygons(face_polygons, polylines):
	pass

def test_dog_from_face_polygons(svg_file):
	print 'Testing with file ', svg_file
	border_poly,polylines = svg_creases_to_polygonal_data(svg_file)
	face_polygons, polylines = crease_pattern(border_poly, polylines)

if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_from_face_polygons(sys.argv[1])
	else:
		test_dog_from_face_polygons("../crease_patterns/1_curve.svg")
