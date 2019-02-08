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
A0[0][4] = -lambda*t30;
  A0[0][5] = lambda*(t3*(t25-ep_0_t*(t11+t15-ep_f_t*ep_f_v1_y-ep_0_v2_y*t3))-ep_0_t*t58);
  A0[0][7] = -lambda*t39;
  A0[0][8] = t67;
  A0[0][10] = t44;
  A0[0][11] = -lambda*t73;
  A0[0][13] = t49;
  A0[0][14] = -lambda*t80;
  A0[0][16] = -lambda*t53;
  A0[0][17] = t88;
  A0[0][19] = -lambda*t24;
  A0[0][20] = t93;
  A0[1][3] = t96;
  A0[1][5] = -lambda*t63;
  A0[1][6] = t149;
  A0[1][8] = -lambda*t70;
  A0[1][9] = -t44;
  A0[1][11] = t77;
  A0[1][12] = -t49;
  A0[1][14] = t84;
  A0[1][15] = t172;
  A0[1][17] = -lambda*t91;
  A0[1][18] = t179;
  A0[1][20] = -lambda*t55;
  A0[2][3] = t99;
  A0[2][4] = t102;
  A0[2][6] = -t67;
  A0[2][7] = t151;
  A0[2][9] = t156;
  A0[2][10] = -t77;
  A0[2][12] = t164;
  A0[2][13] = -t84;
  A0[2][15] = -t88;
  A0[2][16] = t176;
  A0[2][18] = -t93;
  A0[2][19] = t185;
  A0[3][1] = t96;
  A0[3][2] = t99;
  A0[3][7] = t106;
  A0[3][8] = -lambda*t122;
  A0[3][10] = -lambda*t109;
  A0[3][11] = t130;
  A0[3][13] = -lambda*t112;
  A0[3][14] = t137;
  A0[3][16] = t116;
  A0[3][17] = -lambda*t143;
  A0[3][19] = t117;
  A0[3][20] = -lambda*t58;
  A0[4][0] = -t96;
  A0[4][2] = t102;
  A0[4][6] = -t106;
  A0[4][8] = t126;
  A0[4][9] = t157;
  A0[4][11] = -lambda*t133;
  A0[4][12] = t165;
  A0[4][14] = -lambda*t140;
  A0[4][15] = -t116;
  A0[4][17] = t147;
  A0[4][18] = -t117;
  A0[4][20] = t148;
  A0[5][0] = -t99;
  A0[5][1] = -t102;
  A0[5][6] = t150;
  A0[5][7] = -t126;
  A0[5][9] = -t130;
  A0[5][10] = t158;
  A0[5][12] = -t137;
  A0[5][13] = t168;
  A0[5][15] = t173;
  A0[5][16] = -t147;
  A0[5][18] = t180;
  A0[5][19] = -t148;
  A0[6][1] = t149;
  A0[6][2] = -t67;
  A0[6][4] = -t106;
  A0[6][5] = t150;
  A0[6][13] = -ep_b_t*ep_f_t*lambda*t18;
  A0[6][14] = t153;
  A0[6][16] = t152;
  A0[6][17] = -ep_b_t*lambda*t8*t19;
  A0[6][19] = -ep_b_t*lambda*t10;
  A0[6][20] = t155;
  A0[7][0] = -t149;
  A0[7][2] = t151;
  A0[7][3] = t106;
  A0[7][5] = -t126;
  A0[7][12] = t166;
  A0[7][14] = -ep_b_t*ep_f_t*lambda*t45;
  A0[7][15] = -t152;
  A0[7][17] = t154;
  A0[7][18] = t181;
  A0[7][20] = -ep_b_t*lambda*t36;
  A0[8][0] = t67;
  A0[8][1] = -t151;
  A0[8][3] = -t150;
  A0[8][4] = t126;
  A0[8][12] = -t153;
  A0[8][13] = t169;
  A0[8][15] = t174;
  A0[8][16] = -t154;
  A0[8][18] = -t155;
  A0[8][19] = t186;
  A0[9][1] = -t44;
  A0[9][2] = t156;
  A0[9][4] = t157;
  A0[9][5] = -t130;
  A0[9][13] = t159;
  A0[9][14] = -ep_f_t*lambda*t4*t19;
  A0[9][16] = -lambda*t4*t8*t18;
  A0[9][17] = t162;
  A0[9][19] = t160;
  A0[9][20] = -lambda*t4*t17;
  A0[10][0] = t44;
  A0[10][2] = -t77;
  A0[10][3] = -t157;
  A0[10][5] = t158;
  A0[10][12] = -t159;
  A0[10][14] = t161;
  A0[10][15] = t175;
  A0[10][17] = -lambda*t4*t8*t45;
  A0[10][18] = -t160;
  A0[10][20] = t163;
  A0[11][0] = -t156;
  A0[11][1] = t77;
  A0[11][3] = t130;
  A0[11][4] = -t158;
  A0[11][12] = t167;
  A0[11][13] = -t161;
  A0[11][15] = -t162;
  A0[11][16] = t177;
  A0[11][18] = t182;
  A0[11][19] = -t163;
  A0[12][1] = -t49;
  A0[12][2] = t164;
  A0[12][4] = t165;
  A0[12][5] = -t137;
  A0[12][7] = t166;
  A0[12][8] = -t153;
  A0[12][10] = -t159;
  A0[12][11] = t167;
  A0[12][19] = t170;
  A0[12][20] = -ep_f_t*lambda*t13;
  A0[13][0] = t49;
  A0[13][2] = -t84;
  A0[13][3] = -t165;
  A0[13][5] = t168;
  A0[13][6] = -t166;
  A0[13][8] = t169;
  A0[13][9] = t159;
  A0[13][11] = -t161;
  A0[13][18] = -t170;
  A0[13][20] = t171;
  A0[14][0] = -t164;
  A0[14][1] = t84;
  A0[14][3] = t137;
  A0[14][4] = -t168;
  A0[14][6] = t153;
  A0[14][7] = -t169;
  A0[14][9] = -t167;
  A0[14][10] = t161;
  A0[14][18] = t183;
  A0[14][19] = -t171;
  A0[15][1] = t172;
  A0[15][2] = -t88;
  A0[15][4] = -t116;
  A0[15][5] = t173;
  A0[15][7] = -t152;
  A0[15][8] = t174;
  A0[15][10] = t175;
  A0[15][11] = -t162;
  A0[15][19] = -lambda*t6*t8;
  A0[15][20] = t178;
  A0[16][0] = -t172;
  A0[16][2] = t176;
  A0[16][3] = t116;
  A0[16][5] = -t147;
  A0[16][6] = t152;
  A0[16][8] = -t154;
  A0[16][9] = -t175;
  A0[16][11] = t177;
  A0[16][18] = t184;
  A0[16][20] = -lambda*t8*t33;
  A0[17][0] = t88;
  A0[17][1] = -t176;
  A0[17][3] = -t173;
  A0[17][4] = t147;
  A0[17][6] = -t174;
  A0[17][7] = t154;
  A0[17][9] = t162;
  A0[17][10] = -t177;
  A0[17][18] = -t178;
  A0[17][19] = t187;
  A0[18][1] = t179;
  A0[18][2] = -t93;
  A0[18][4] = -t117;
  A0[18][5] = t180;
  A0[18][7] = t181;
  A0[18][8] = -t155;
  A0[18][10] = -t160;
  A0[18][11] = t182;
  A0[18][13] = -t170;
  A0[18][14] = t183;
  A0[18][16] = t184;
  A0[18][17] = -t178;
  A0[19][0] = -t179;
  A0[19][2] = t185;
  A0[19][3] = t117;
  A0[19][5] = -t148;
  A0[19][6] = -t181;
  A0[19][8] = t186;
  A0[19][9] = t160;
  A0[19][11] = -t163;
  A0[19][12] = t170;
  A0[19][14] = -t171;
  A0[19][15] = -t184;
  A0[19][17] = t187;
  A0[20][0] = t93;
  A0[20][1] = -t185;
  A0[20][3] = -t180;
  A0[20][4] = t148;
  A0[20][6] = t155;
  A0[20][7] = -t186;
  A0[20][9] = -t182;
  A0[20][10] = t163;
  A0[20][12] = -t183;
  A0[20][13] = t171;
  A0[20][15] = t178;
  A0[20][16] = -t187;
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
print matlab_hess_to_c_code(MV_const_dict, hessian_str)
#print double_me(double_me_str)
#print matlab_cross_grad_to_c_code(grad_str)
#print matlab_cross_hessian_to_c_code(hessian_str)
#print matlab_cross_bnd_hessian_to_c_code(hessian_str)
#print matlab_edge_hessian_to_c_code(hessian_str)