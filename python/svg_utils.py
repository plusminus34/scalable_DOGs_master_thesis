from svgpathtools import svg2paths, wsvg
from xml.dom import minidom
import struct
import numpy as np

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
			print 'cls = ', cls
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
	return struct.unpack('BBB',color_undecoded.strip('#').decode('hex'))

def getText(nodelist):                                              
	rc = []
	for node in nodelist:
		if node.nodeType == node.TEXT_NODE:
			rc.append(node.data)
	return ''.join(rc)


def get_path_points_as_matrix(path, points_num, center_x, center_y):
	points = np.array([path.point(t) for t in np.arange(0,1,1./points_num)])
	points = points.view('(2,)float')
	points[:,1] = -1*(points[:,1]-center_y)
	points[:,1] = points[:,1] - center_y
	return points