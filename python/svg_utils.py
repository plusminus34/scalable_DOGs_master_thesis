from svgpathtools import svg2paths, wsvg
from xml.dom import minidom
import struct
import numpy as np
import binascii
from shapely.ops import polygonize
"""
>>> from shapely.ops import polygonize
>>> lines = [
...     ((0, 0), (1, 1)),
...     ((0, 0), (0, 1)),
...     ((0, 1), (1, 1)),
...     ((1, 1), (1, 0)),
...     ((1, 0), (0, 0))
...     ]
"""

def get_style_classes(svg_file):
	doc = minidom.parse(svg_file)
	if (len(doc.getElementsByTagName("style")) < 1):
		return {}
	txt = getText(doc.getElementsByTagName("style")[0].childNodes).split('\n')
	style_classes = {}
	for i in range(1,len(txt)-1):
		style_cls = {}
		txt_l = txt[i]
		#print 'txt_l = ', txt_l
		txt_f = txt_l[txt_l.find('.')+1:]
		#print 'txt_f = ', txt_f
		cls = txt_f[0:txt_f.find('{')]
		#print 'cls = ', cls
		content = txt_f[txt_f.find('{')+1:txt_f.find('}')-1].split(';')
		#print 'content = ', content
		for cont in content:
			cont_splt = cont.split(":")
			#print 'cont_splt = ',cont_splt
			style_cls[cont_splt[0]]= cont_splt[1]
		#print 'style_cls = ', style_cls
		style_classes[cls] = style_cls

	#print 'txt = ', txt
	return style_classes

def get_svg_view_box(svg_file):
	doc = minidom.parse(svg_file)
	viewBoxStr = doc.getElementsByTagName("svg")[0].getAttribute('viewBox')
	viewBoxList = viewBoxStr.split(' ')
	viewBoxList[1] = str(float(viewBoxList[1]) - float(viewBoxList[3]))
	return [float(x) for x in viewBoxList]

def get_curve_color(style_classes,attrib):
	#print 'attrib = ', attrib
	#print "attrib[class'] = ", attrib['class']
	if (style_classes != {}):
		# There's a class attribute
		#print "attrib[class'] = ", attrib['class']
		try:
			cls = attrib['class']
			#print 'cls = ', cls
			try:
				color_undecoded = style_classes[cls]['stroke']
			except:
				color_undecoded = style_classes[cls]['fill']
			#print 'color_undecoded = ', color_undecoded
		except:
			color_undecoded = "#000000"
	else:
		color_undecoded = attrib['stroke']
	# The color information is in the path
	return struct.unpack( 'BBB', binascii.unhexlify(color_undecoded.strip('#')) )

def getText(nodelist):                                              
	rc = []
	for node in nodelist:
		if node.nodeType == node.TEXT_NODE:
			rc.append(node.data)
	return ''.join(rc)

def sample_bezier_path_sampling(path, points_num):
	points = np.array([path.point(t) for t in np.arange(0,1+1./points_num,1./points_num)])
	points = points.view('(2,)float')
	#points[:,1] = -1*(points[:,1]-center_y)
	#points[:,1] = points[:,1] - center_y
	#print 'points.shape = ', points.shape
	return points

def sample_polylines(path,bounds):
	#print 'len(path) = ', len(path)
	print ('sampling a polyline!')
	points = np.empty((len(path)+1,2))
	#print 'path[0].start.real,path[0].start.imag = ', path[0].start.real,path[0].start.imag
	points[0,:] = path[0].start.real,path[0].start.imag
	idx = 1
	for line in path:
		points[idx,:] = line.end.real,line.end.imag
		#print 'points[idx,:] = ', points[idx,:]
		idx += 1
	#print 'points.shape = ', points.shape
	#print 'points = ', points
	
	eps = 1e-2*abs(bounds[2]-bounds[0])
	points[0,0] = round_if_close(points[0,0], bounds[0],eps)
	points[0,0] = round_if_close(points[0,0], bounds[2],eps)
	points[0,1] = round_if_close(points[0,1], -bounds[1],eps)
	points[0,1] = round_if_close(points[0,1], -bounds[3],eps)

	points[-1,0] = round_if_close(points[-1,0], bounds[0],eps)
	points[-1,0] = round_if_close(points[-1,0], bounds[2],eps)
	points[-1,1] = round_if_close(points[-1,1],-bounds[1],eps)
	points[-1,1] = round_if_close(points[-1,1],-bounds[3],eps)

	return points

def round_if_close(pt,close_pt,eps):
	print ('pt = ', pt, ' close_pt = ', close_pt)
	if (abs(pt-close_pt) < eps):
		print ('rounding to ', close_pt)
		return close_pt
	return pt

def is_border(attrib,style_classes):
	color = get_curve_color(style_classes,attrib)
	#print 'The color is ', color
	is_border = (color == (0,0,0))
	return is_border
	
	return points
