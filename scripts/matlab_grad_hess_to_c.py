import re
import itertools

fold_bias_var_dict = {}
fold_bias_var_dict['[0]'] = 'ep_b_v1_i'
fold_bias_var_dict['[1]'] = 'ep_b_v1_i+vnum'
fold_bias_var_dict['[2]'] = 'ep_b_v1_i+2*vnum'
fold_bias_var_dict['[3]'] = 'ep_b_v2_i'
fold_bias_var_dict['[4]'] = 'ep_b_v2_i+vnum'
fold_bias_var_dict['[5]'] = 'ep_b_v2_i+2*vnum'
fold_bias_var_dict['[6]'] = 'ep_f_v1_i'
fold_bias_var_dict['[7]'] = 'ep_f_v1_i+vnum'
fold_bias_var_dict['[8]'] = 'ep_f_v1_i+2*vnum'
fold_bias_var_dict['[9]'] = 'ep_f_v2_i'
fold_bias_var_dict['[10]'] = 'ep_f_v2_i+vnum'
fold_bias_var_dict['[11]'] = 'ep_f_v2_i+2*vnum'
fold_bias_var_dict['[12]'] = 'v1_i'
fold_bias_var_dict['[13]'] = 'v1_i+vnum'
fold_bias_var_dict['[14]'] = 'v1_i+2*vnum'
fold_bias_var_dict['[15]'] = 'v2_i'
fold_bias_var_dict['[16]'] = 'v2_i+vnum'
fold_bias_var_dict['[17]'] = 'v2_i+2*vnum'
fold_bias_var_dict['[18]'] = 'w1_i'
fold_bias_var_dict['[19]'] = 'w1_i+vnum'
fold_bias_var_dict['[20]'] = 'w1_i+2*vnum'
fold_bias_var_dict['[21]'] = 'w2_i'
fold_bias_var_dict['[22]'] = 'w2_i+vnum'
fold_bias_var_dict['[23]'] = 'w2_i+2*vnum'



grad_str = """
  A0[0][0] = t48*(ep_b_t*t6*t10-ep_b_t*t9*t15+ep_b_t*t6*t26-ep_b_t*t9*t31)*-2.0;
  A0[1][0] = t48*(ep_b_t*t10*t30-ep_b_t*t9*t39-ep_b_t*t9*t40+ep_b_t*t26*t30)*2.0;
  A0[2][0] = t48*(ep_b_t*t6*t39-ep_b_t*t15*t30+ep_b_t*t6*t40-ep_b_t*t30*t31)*2.0;
  A0[3][0] = t48*(t6*t10*t13-t9*t13*t15+t6*t13*t26-t9*t13*t31)*2.0;
  A0[4][0] = t48*(t10*t13*t30-t9*t13*t39-t9*t13*t40+t13*t26*t30)*-2.0;
  A0[5][0] = t48*(t6*t13*t39-t13*t15*t30+t6*t13*t40-t13*t30*t31)*-2.0;
  A0[6][0] = (ep_f_t*t10*t23-ep_f_t*t15*t25+ep_f_t*t23*t26-ep_f_t*t25*t31)*(t43+t45+t46+t47-t49-t50)*2.0;
  A0[7][0] = t48*(ep_f_t*t10*t17+ep_f_t*t17*t26-ep_f_t*t25*t39-ep_f_t*t25*t40)*-2.0;
  A0[8][0] = (ep_f_t*t15*t17+ep_f_t*t17*t31-ep_f_t*t23*t39-ep_f_t*t23*t40)*(t43+t45+t46+t47-t49-t50)*2.0;
  A0[9][0] = t48*(t2*t10*t23-t2*t15*t25+t2*t23*t26-t2*t25*t31)*-2.0;
  A0[10][0] = (t2*t10*t17+t2*t17*t26-t2*t25*t39-t2*t25*t40)*(t43+t45+t46+t47-t49-t50)*2.0;
  A0[11][0] = t48*(t2*t15*t17+t2*t17*t31-t2*t23*t39-t2*t23*t40)*-2.0;
  A0[12][0] = (t43+t45+t46+t47-t49-t50)*(t36-t44+t10*t52+t26*t52-t15*(t53-ep_0_t*t25)-t31*(t53-ep_0_t*t25))*2.0;
  A0[13][0] = (t43+t45+t46+t47-t49-t50)*(t32-t42+t10*t55+t26*t55+t39*(t53-ep_0_t*t25)+t40*(t53-ep_0_t*t25))*2.0;
  A0[14][0] = t48*(t28-t41+t15*t55+t31*t55+t39*t52+t40*t52)*-2.0;
  A0[15][0] = t48*(t36-t44+t10*t59-t15*t61+t26*t59-t31*t61)*-2.0;
  A0[16][0] = t48*(t32-t42+t10*t64+t26*t64+t39*t61+t40*t61)*-2.0;
  A0[17][0] = (t28-t41+t15*t64+t31*t64+t39*(t58-t66)+t40*(t58-t66))*(t43+t45+t46+t47-t49-t50)*2.0;
  A0[18][0] = t71;
  A0[19][0] = t35*t70*-2.0;
  A0[20][0] = t72;
  A0[21][0] = -t71;
  A0[22][0] = t35*t70*2.0;
  A0[23][0] = -t72;
"""

def matlab_grad_to_c_code(var_dict, str):
	for var_idx in var_dict:
		str = str.replace('A0' + var_idx + '[0] ', 'grad(' + var_dict[var_idx] + ') +')
	return str


def matlab_hess_to_c_code(var_dict, str):
	lines = str.split('\n')
	#print 'lines = ', lines
	code = ''
	#print 'var_dict.values() = ',var_dict.values()
	i = 0
	#order_prod = list(itertools.product(var_orders, var_orders))
	for l in lines:
		if len(l) > 10:
			first_part = l[l.find('['):l.find(']')+1]
			second_start_idx = l.find(']')+1
			l_sec = l[second_start_idx:]
			second_part = l_sec[0:l_sec.find(']')+1]
			g1_i,g2_i = var_dict[first_part],var_dict[second_part]
			val = ''.join(l.split('=')[1:]).strip(';')
			#print 'val = ', val
			#sys.exit(1)
			#g1_i,g2_i = order_prod[i][0],order_prod[i][1]
			new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			code += new_t
			i+=1
	return code

hessian_str = """
A0[12][12] = t3;
  A0[12][13] = t4;
  A0[12][14] = t5;
  A0[12][15] = -t3;
  A0[12][16] = -t4;
  A0[12][17] = -t5;
  A0[12][18] = t3;
  A0[12][19] = t4;
  A0[12][20] = t5;
  A0[12][21] = -t3;
  A0[12][22] = -t4;
  A0[12][23] = -t5;
  A0[13][12] = t4;
  A0[13][13] = t7;
  A0[13][14] = t8;
  A0[13][15] = -t4;
  A0[13][16] = -t7;
  A0[13][17] = -t8;
  A0[13][18] = t4;
  A0[13][19] = t7;
  A0[13][20] = t8;
  A0[13][21] = -t4;
  A0[13][22] = -t7;
  A0[13][23] = -t8;
  A0[14][12] = t5;
  A0[14][13] = t8;
  A0[14][14] = t10;
  A0[14][15] = -t5;
  A0[14][16] = -t8;
  A0[14][17] = -t10;
  A0[14][18] = t5;
  A0[14][19] = t8;
  A0[14][20] = t10;
  A0[14][21] = -t5;
  A0[14][22] = -t8;
  A0[14][23] = -t10;
  A0[15][12] = -t3;
  A0[15][13] = -t4;
  A0[15][14] = -t5;
  A0[15][15] = t3;
  A0[15][16] = t4;
  A0[15][17] = t5;
  A0[15][18] = -t3;
  A0[15][19] = -t4;
  A0[15][20] = -t5;
  A0[15][21] = t3;
  A0[15][22] = t4;
  A0[15][23] = t5;
  A0[16][12] = -t4;
  A0[16][13] = -t7;
  A0[16][14] = -t8;
  A0[16][15] = t4;
  A0[16][16] = t7;
  A0[16][17] = t8;
  A0[16][18] = -t4;
  A0[16][19] = -t7;
  A0[16][20] = -t8;
  A0[16][21] = t4;
  A0[16][22] = t7;
  A0[16][23] = t8;
  A0[17][12] = -t5;
  A0[17][13] = -t8;
  A0[17][14] = -t10;
  A0[17][15] = t5;
  A0[17][16] = t8;
  A0[17][17] = t10;
  A0[17][18] = -t5;
  A0[17][19] = -t8;
  A0[17][20] = -t10;
  A0[17][21] = t5;
  A0[17][22] = t8;
  A0[17][23] = t10;
  A0[18][12] = t3;
  A0[18][13] = t4;
  A0[18][14] = t5;
  A0[18][15] = -t3;
  A0[18][16] = -t4;
  A0[18][17] = -t5;
  A0[18][18] = t3;
  A0[18][19] = t4;
  A0[18][20] = t5;
  A0[18][21] = -t3;
  A0[18][22] = -t4;
  A0[18][23] = -t5;
  A0[19][12] = t4;
  A0[19][13] = t7;
  A0[19][14] = t8;
  A0[19][15] = -t4;
  A0[19][16] = -t7;
  A0[19][17] = -t8;
  A0[19][18] = t4;
  A0[19][19] = t7;
  A0[19][20] = t8;
  A0[19][21] = -t4;
  A0[19][22] = -t7;
  A0[19][23] = -t8;
  A0[20][12] = t5;
  A0[20][13] = t8;
  A0[20][14] = t10;
  A0[20][15] = -t5;
  A0[20][16] = -t8;
  A0[20][17] = -t10;
  A0[20][18] = t5;
  A0[20][19] = t8;
  A0[20][20] = t10;
  A0[20][21] = -t5;
  A0[20][22] = -t8;
  A0[20][23] = -t10;
  A0[21][12] = -t3;
  A0[21][13] = -t4;
  A0[21][14] = -t5;
  A0[21][15] = t3;
  A0[21][16] = t4;
  A0[21][17] = t5;
  A0[21][18] = -t3;
  A0[21][19] = -t4;
  A0[21][20] = -t5;
  A0[21][21] = t3;
  A0[21][22] = t4;
  A0[21][23] = t5;
  A0[22][12] = -t4;
  A0[22][13] = -t7;
  A0[22][14] = -t8;
  A0[22][15] = t4;
  A0[22][16] = t7;
  A0[22][17] = t8;
  A0[22][18] = -t4;
  A0[22][19] = -t7;
  A0[22][20] = -t8;
  A0[22][21] = t4;
  A0[22][22] = t7;
  A0[22][23] = t8;
  A0[23][12] = -t5;
  A0[23][13] = -t8;
  A0[23][14] = -t10;
  A0[23][15] = t5;
  A0[23][16] = t8;
  A0[23][17] = t10;
  A0[23][18] = -t5;
  A0[23][19] = -t8;
  A0[23][20] = -t10;
  A0[23][21] = t5;
  A0[23][22] = t8;
  A0[23][23] = t10;
 """

def matlab_cross_hessian_to_c_code(str):
	var_dict = {}
	var_dict['[0]'] = 'p_0_i'
	var_dict['[1]'] = 'p_0_i+vnum'
	var_dict['[2]'] = 'p_0_i+2*vnum'
	var_dict['[3]'] = 'p_xb_i'
	var_dict['[4]'] = 'p_xb_i+vnum'
	var_dict['[5]'] = 'p_xb_i+2*vnum'
	var_dict['[6]'] = 'p_xf_i'
	var_dict['[7]'] = 'p_xf_i+vnum'
	var_dict['[8]'] = 'p_xf_i+2*vnum'
	var_dict['[9]'] = 'p_yb_i'
	var_dict['[10]'] = 'p_yb_i+vnum'
	var_dict['[11]'] = 'p_yb_i+2*vnum'
	var_dict['[12]'] = 'p_yf_i'
	var_dict['[13]'] = 'p_yf_i+vnum'
	var_dict['[14]'] = 'p_yf_i+2*vnum'

	#var_orders = ['p_0_i','p_0_i+vnum','p_0_i+2*vnum','p_xb_i','p_xb_i+vnum','p_xb_i+2*vnum', \
	#				'p_xf_i', 'p_xf_i+vnum','p_xf_i+2*vnum','p_yb_i','p_yb_i+vnum','p_yb_i+2*vnum', \
	#				'p_yf_i','p_yf_i+vnum','p_yf_i+2*vnum']
	#order_prod = list(itertools.product(var_orders, var_orders))

	lines = str.split('\n')
	#print 'lines = ', lines
	code = ''
	#print 'var_dict.values() = ',var_dict.values()
	i = 0
	#order_prod = list(itertools.product(var_orders, var_orders))
	for l in lines:
		if len(l) > 10:
			first_part = l[l.find('['):l.find(']')+1]
			second_start_idx = l.find(']')+1
			l_sec = l[second_start_idx:]
			second_part = l_sec[0:l_sec.find(']')+1]
			g1_i,g2_i = var_dict[first_part],var_dict[second_part]
			val = ''.join(l.split('=')[1:]).strip(';')
			#print 'val = ', val
			#sys.exit(1)
			#g1_i,g2_i = order_prod[i][0],order_prod[i][1]
			new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			code += new_t
			i+=1
	return code

def matlab_cross_bnd_hessian_to_c_code(str):
	var_dict = {}
	var_dict['[0]'] = 'p_0_i'
	var_dict['[1]'] = 'p_0_i+vnum'
	var_dict['[2]'] = 'p_0_i+2*vnum'
	var_dict['[3]'] = 'p_xb_i'
	var_dict['[4]'] = 'p_xb_i+vnum'
	var_dict['[5]'] = 'p_xb_i+2*vnum'
	var_dict['[6]'] = 'p_xf_i'
	var_dict['[7]'] = 'p_xf_i+vnum'
	var_dict['[8]'] = 'p_xf_i+2*vnum'
	var_dict['[9]'] = 'p_yf_i'
	var_dict['[10]'] = 'p_yf_i+vnum'
	var_dict['[11]'] = 'p_yf_i+2*vnum'

	#var_orders = ['p_0_i','p_0_i+vnum','p_0_i+2*vnum','p_xb_i','p_xb_i+vnum','p_xb_i+2*vnum', \
	#				'p_xf_i', 'p_xf_i+vnum','p_xf_i+2*vnum', \
  #				'p_yf_i','p_yf_i+vnum','p_yf_i+2*vnum']
	#order_prod = list(itertools.product(var_orders, var_orders))

	lines = str.split('\n')
	#print 'lines = ', lines
	code = ''
	#print 'var_dict.values() = ',var_dict.values()
	i = 0
	#order_prod = list(itertools.product(var_orders, var_orders))
	for l in lines:
		if len(l) > 10:
			first_part = l[l.find('['):l.find(']')+1]
			second_start_idx = l.find(']')+1
			l_sec = l[second_start_idx:]
			second_part = l_sec[0:l_sec.find(']')+1]
			g1_i,g2_i = var_dict[first_part],var_dict[second_part]
			val = ''.join(l.split('=')[1:]).strip(';')
			new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			code += new_t
			i+=1
	return code

#[ p0_x, p0_y, p0_z, pxb_x, pxb_y, pxb_z, pxf_x, pxf_y, pxf_z]
def matlab_curve_hessian_to_c_code(str):
	var_dict = {}
	var_dict['[0]'] = '[p_0_i]'
	var_dict['[1]'] = '[p_0_i+vnum]'
	var_dict['[2]'] = '[p_0_i+2*vnum]'
	var_dict['[3]'] = '[p_xb_i]'
	var_dict['[4]'] = '[p_xb_i+vnum]'
	var_dict['[5]'] = '[p_xb_i+2*vnum]'
	var_dict['[6]'] = '[p_xf_i]'
	var_dict['[7]'] = '[p_xf_i+vnum]'
	var_dict['[8]'] = '[p_xf_i+2*vnum]'

	var_orders = ['p_0_i','p_0_i+vnum','p_0_i+2*vnum','p_xb_i','p_xb_i+vnum','p_xb_i+2*vnum', \
					'p_xf_i', 'p_xf_i+vnum','p_xf_i+2*vnum']
	order_prod = list(itertools.product(var_orders, var_orders))

	lines = str.split('\n')
	code = ''
	#print 'var_dict.values() = ',var_dict.values()
	i = 0
	order_prod = list(itertools.product(var_orders, var_orders))
	for l in lines:
		if len(l) > 10:
			val = ''.join(l.split('=')[1:]).strip(';')
			g1_i,g2_i = order_prod[i][0],order_prod[i][1]
			new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			code += new_t
			i+=1

	#print 'processed ', i , ' lines'
	return code


double_me_str = """  t2 = p0_x-pxf_x;
	t3 = p0_y-pxf_y;
	t4 = p0_z-pxf_z;
	t5 = p0_x-pyb_x;
	t6 = p0_y-pyb_y;
	t7 = p0_z-pyb_z;
	t8 = p0_x-pxb_x;
	t9 = p0_y-pxb_y;
	t10 = p0_z-pxb_z;
	t11 = t5*t8;
	t12 = t6*t9;
	t97 = t85-t96;
	t98 = pyf_y*2.0;
	t99 = t85-t98;
	t100 = t22*t62;
	t575 = t91*(t384-t399)*2.0;
	t932 = t190-t358;"""
def double_me(s):
	lines = s.split('\n')
	#print 'lines = ', lines
	code = ''
	#print 'var_dict.values() = ',var_dict.values()
	i = 2 # matlab's output starts with 2 for some reason
	for l in lines:
		if len(l) > 2:
			#print 'l = ', l
			val = ''.join(l.split('=')[1:]).strip(';')
			#print 'val = ', val
			new_t = 'double t' + str(int(i)) + ' = ' + val + ';\n'
			#print 'new_t = ', new_t
			#new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			if new_t:
				code += new_t
				i+=1
	return code

#print matlab_grad_to_c_code(fold_bias_var_dict, grad_str)
print matlab_hess_to_c_code(fold_bias_var_dict, hessian_str)
#print double_me(double_me_str)
#print matlab_cross_grad_to_c_code(grad_str)
#print matlab_cross_hessian_to_c_code(hessian_str)
#print matlab_cross_bnd_hessian_to_c_code(hessian_str)
#print matlab_edge_hessian_to_c_code(hessian_str)