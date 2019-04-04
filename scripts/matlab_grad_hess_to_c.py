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

dog_star_var_dict = {}
dog_star_var_dict['[0]'] = 'p_0_i'
dog_star_var_dict['[1]'] = 'p_0_i+vnum'
dog_star_var_dict['[2]'] = 'p_0_i+2*vnum'
dog_star_var_dict['[3]'] = 'p_xb_i'
dog_star_var_dict['[4]'] = 'p_xb_i+vnum'
dog_star_var_dict['[5]'] = 'p_xb_i+2*vnum'
dog_star_var_dict['[6]'] = 'p_xf_i'
dog_star_var_dict['[7]'] = 'p_xf_i+vnum'
dog_star_var_dict['[8]'] = 'p_xf_i+2*vnum'
dog_star_var_dict['[9]'] = 'p_yb_i'
dog_star_var_dict['[10]'] = 'p_yb_i+vnum'
dog_star_var_dict['[11]'] = 'p_yb_i+2*vnum'
dog_star_var_dict['[12]'] = 'p_yf_i'
dog_star_var_dict['[13]'] = 'p_yf_i+vnum'
dog_star_var_dict['[14]'] = 'p_yf_i+2*vnum'
 
dog_star_var_dict_bnd = {}
dog_star_var_dict_bnd['[0]'] = 'p_0_i'
dog_star_var_dict_bnd['[1]'] = 'p_0_i+vnum'
dog_star_var_dict_bnd['[2]'] = 'p_0_i+2*vnum'
dog_star_var_dict_bnd['[3]'] = 'p_xb_i'
dog_star_var_dict_bnd['[4]'] = 'p_xb_i+vnum'
dog_star_var_dict_bnd['[5]'] = 'p_xb_i+2*vnum'
dog_star_var_dict_bnd['[6]'] = 'p_xf_i'
dog_star_var_dict_bnd['[7]'] = 'p_xf_i+vnum'
dog_star_var_dict_bnd['[8]'] = 'p_xf_i+2*vnum'
dog_star_var_dict_bnd['[9]'] = 'p_yf_i'
dog_star_var_dict_bnd['[10]'] = 'p_yf_i+vnum'
dog_star_var_dict_bnd['[11]'] = 'p_yf_i+2*vnum'

MV_const_dict = {}
MV_const_dict['[0]'] = 'ep_0_v1_i'
MV_const_dict['[1]'] = 'ep_0_v1_i+vnum'
MV_const_dict['[2]'] = 'ep_0_v1_i+2*vnum'
MV_const_dict['[3]'] = 'ep_0_v2_i'
MV_const_dict['[4]'] = 'ep_0_v2_i+vnum'
MV_const_dict['[5]'] = 'ep_0_v2_i+2*vnum'
MV_const_dict['[6]'] = 'ep_b_v1_i'
MV_const_dict['[7]'] = 'ep_b_v1_i+vnum'
MV_const_dict['[8]'] = 'ep_b_v1_i+2*vnum'
MV_const_dict['[9]'] = 'ep_b_v2_i'
MV_const_dict['[10]'] = 'ep_b_v2_i+vnum'
MV_const_dict['[11]'] = 'ep_b_v2_i+2*vnum'
MV_const_dict['[12]'] = 'ep_f_v1_i'
MV_const_dict['[13]'] = 'ep_f_v1_i+vnum'
MV_const_dict['[14]'] = 'ep_f_v1_i+2*vnum'
MV_const_dict['[15]'] = 'ep_f_v2_i'
MV_const_dict['[16]'] = 'ep_f_v2_i+vnum'
MV_const_dict['[17]'] = 'ep_f_v2_i+2*vnum'
MV_const_dict['[18]'] = 'folded_v_i'
MV_const_dict['[19]'] = 'folded_v_i+vnum'
MV_const_dict['[20]'] = 'folded_v_i+2*vnum'

pair_rigid_align_dict = {}
pair_rigid_align_dict['[0]'] = 'p_1_i'
pair_rigid_align_dict['[1]'] = 'p_1_i+vnum'
pair_rigid_align_dict['[2]'] = 'p_1_i+2*vnum'
pair_rigid_align_dict['[3]'] = 'p_2_i'
pair_rigid_align_dict['[4]'] = 'p_2_i+vnum'
pair_rigid_align_dict['[5]'] = 'p_2_i+2*vnum'

grad_str = """
	A0[0][0] = t56*(alpha*t58*(ep_b_t*t15*t36-ep_b_t*t22*t35)+alpha*t60*(ep_b_t*t15*t47-ep_b_t*t22*t46))*2.0;
	A0[1][0] = t56*(alpha*t58*(ep_b_t*t15*t40-ep_b_t*t20*t35)+alpha*t60*(ep_b_t*t15*t49-ep_b_t*t20*t46))*-2.0;
	A0[2][0] = t56*(alpha*t58*(ep_b_t*t20*t36-ep_b_t*t22*t40)+alpha*t60*(ep_b_t*t20*t47-ep_b_t*t22*t49))*-2.0;
	A0[3][0] = t56*(alpha*t58*(t3*t15*t36-t3*t22*t35)+alpha*t60*(t3*t15*t47-t3*t22*t46))*-2.0;
	A0[4][0] = t56*(alpha*t58*(t3*t15*t40-t3*t20*t35)+alpha*t60*(t3*t15*t49-t3*t20*t46))*2.0;
	A0[5][0] = t56*(alpha*t58*(t3*t20*t36-t3*t22*t40)+alpha*t60*(t3*t20*t47-t3*t22*t49))*2.0;
	A0[6][0] = t56*(alpha*t58*(ep_f_t*t13*t35-ep_f_t*t21*t36)+alpha*t60*(ep_f_t*t13*t46-ep_f_t*t21*t47))*2.0;
	A0[7][0] = t56*(alpha*t58*(ep_f_t*t8*t35-ep_f_t*t21*t40)+alpha*t60*(ep_f_t*t8*t46-ep_f_t*t21*t49))*-2.0;
	A0[8][0] = t56*(alpha*t58*(ep_f_t*t8*t36-ep_f_t*t13*t40)+alpha*t60*(ep_f_t*t8*t47-ep_f_t*t13*t49))*2.0;
	A0[9][0] = t56*(alpha*t58*(t5*t13*t35-t5*t21*t36)+alpha*t60*(t5*t13*t46-t5*t21*t47))*-2.0;
	A0[10][0] = t56*(alpha*t58*(t5*t8*t35-t5*t21*t40)+alpha*t60*(t5*t8*t46-t5*t21*t49))*2.0;
	A0[11][0] = t56*(alpha*t58*(t5*t8*t36-t5*t13*t40)+alpha*t60*(t5*t8*t47-t5*t13*t49))*-2.0;
	A0[12][0] = t56*(alpha*t60*(t46*t62+t47*t64)+alpha*t58*(t33-t41+t35*t62+t36*t64))*-2.0;
	A0[13][0] = t56*(alpha*t60*(t46*t67+t49*(t63-t65))+alpha*t58*(t30-t39+t35*t67+t40*(t63-t65)))*2.0;
	A0[14][0] = t56*(alpha*t60*(t49*t62-t47*t67)-alpha*t58*(t25-t37-t40*t62+t36*t67))*2.0;
	A0[15][0] = t56*(alpha*t60*(t46*t71+t47*t75)+alpha*t58*(t33-t41+t35*t71+t36*t75))*2.0;
	A0[16][0] = t56*(alpha*t60*(t46*t74+t49*t75)+alpha*t58*(t30-t39+t35*t74+t40*t75))*-2.0;
	A0[17][0] = t56*(alpha*t60*(t49*t71-t47*t74)-alpha*t58*(t25-t37+t36*t74-t40*t71))*-2.0;
	A0[18][0] = alpha*t34*t56*t60*-2.0;
	A0[19][0] = t79;
	A0[20][0] = alpha*t28*t56*t60*-2.0;
	A0[21][0] = alpha*t34*t56*t60*2.0;
	A0[22][0] = -t79;
	A0[23][0] = alpha*t28*t56*t60*2.0;
"""

def matlab_grad_to_c_code(var_dict, str):
	for var_idx in var_dict:
		str = str.replace('A0' + var_idx + '[0] ', 'grad(' + var_dict[var_idx] + ') +')
	return str


def matlab_hess_to_c_code(var_dict, str, cached=True):
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
			if cached:
			 new_t = 'IJV[ijv_idx++] = Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+');\n' 
			else:
			 new_t = 'IJV.push_back(Eigen::Triplet<double>('+g1_i+','+g2_i+','+val+'));\n'
			code += new_t
			i+=1
	return code

hessian_str = """
  A0[0][0] = (r11*r11)*2.0+(r12*r12)*2.0+(r13*r13)*2.0;
  A0[0][1] = t5;
  A0[0][2] = t9;
  A0[0][3] = r11*-2.0;
  A0[0][4] = r12*-2.0;
  A0[0][5] = r13*-2.0;
  A0[1][0] = t5;
  A0[1][1] = (r21*r21)*2.0+(r22*r22)*2.0+(r23*r23)*2.0;
  A0[1][2] = t13;
  A0[1][3] = r21*-2.0;
  A0[1][4] = r22*-2.0;
  A0[1][5] = r23*-2.0;
  A0[2][0] = t9;
  A0[2][1] = t13;
  A0[2][2] = (r31*r31)*2.0+(r32*r32)*2.0+(r33*r33)*2.0;
  A0[2][3] = r31*-2.0;
  A0[2][4] = r32*-2.0;
  A0[2][5] = r33*-2.0;
  A0[3][0] = r11*-2.0;
  A0[3][1] = r21*-2.0;
  A0[3][2] = r31*-2.0;
  A0[3][3] = 2.0;
  A0[4][0] = r12*-2.0;
  A0[4][1] = r22*-2.0;
  A0[4][2] = r32*-2.0;
  A0[4][4] = 2.0;
  A0[5][0] = r13*-2.0;
  A0[5][1] = r23*-2.0;
  A0[5][2] = r33*-2.0;
  A0[5][5] = 2.0;
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
print matlab_hess_to_c_code(pair_rigid_align_dict, hessian_str)
#print double_me(double_me_str)
#print matlab_cross_grad_to_c_code(grad_str)
#print matlab_cross_hessian_to_c_code(hessian_str)
#print matlab_cross_bnd_hessian_to_c_code(hessian_str)
#print matlab_edge_hessian_to_c_code(hessian_str)