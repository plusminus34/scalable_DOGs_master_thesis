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
A0[0][0] = -lambda1*(t70+t77+t80+t95-t15*t29*2.0-t20*t24*t30-t15*t24*t43+t20*t30*t40*t46+t15*t35*t43*t46-t29*t30*t40*t47-t15*t45*t47*t48-t29*t34*t41*t81*(3.0/4.0)+t20*t24*(t38*t38)*t41*(3.0/4.0)-t15*t34*(t37*t37)*t59*(3.0/4.0)+t24*t30*t35*t40*t43*(1.0/2.0)-t30*t34*t40*t45*t48*(1.0/2.0));
  A0[0][1] = -lambda1*(t71+t72+t73+t84+t85+t86+t90+t92-t29*t30*t40*t53*(1.0/2.0)-t29*t30*t47*t52*(1.0/2.0)-t15*t45*t48*t53*(1.0/2.0)-t15*t45*t47*t58*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0)-t30*t34*t40*t45*t58*(1.0/4.0)-t30*t34*t45*t48*t52*(1.0/4.0)-t15*t34*t48*t58*t59*(3.0/4.0));
  A0[0][2] = -lambda1*(t74+t75+t76+t87+t88+t89+t109+t111-t29*t30*t40*t64*(1.0/2.0)-t29*t30*t47*t63*(1.0/2.0)-t15*t45*t48*t64*(1.0/2.0)-t15*t45*t47*t68*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0)-t30*t34*t40*t45*t68*(1.0/4.0)-t30*t34*t45*t48*t63*(1.0/4.0)-t15*t34*t48*t59*t68*(3.0/4.0));
  A0[0][3] = t128;
  A0[0][4] = lambda1*(t71+t72+t73-t6*t20*t30*t40*(1.0/2.0)-t6*t15*t35*t43*(1.0/2.0));
  A0[0][5] = lambda1*(t74+t75+t76-t7*t20*t30*t40*(1.0/2.0)-t7*t15*t35*t43*(1.0/2.0));
  A0[0][6] = -lambda1*t157;
  A0[0][7] = -lambda1*(t91+t93+t94-t6*t29*t30*t40*(1.0/2.0)-t6*t15*t45*t48*(1.0/2.0));
  A0[0][8] = -lambda1*(t110+t112+t113-t7*t29*t30*t40*(1.0/2.0)-t7*t15*t45*t48*(1.0/2.0));
  A0[0][9] = t201;
  A0[0][10] = lambda1*(t84+t85+t86+t168+t221-t3*t20*t30*t40*(1.0/2.0)-t3*t15*t35*t43*(1.0/2.0)-t29*t30*t47*t52*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0)-t30*t34*t45*t48*t52*(1.0/4.0));
  A0[0][11] = lambda1*(t87+t88+t89+t170+t222-t4*t20*t30*t40*(1.0/2.0)-t4*t15*t35*t43*(1.0/2.0)-t29*t30*t47*t63*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0)-t30*t34*t45*t48*t63*(1.0/4.0));
  A0[1][0] = -lambda1*(t71+t72+t73+t84+t85+t86+t90-t91+t92-t93-t94-t29*t30*t40*t53*(1.0/2.0)-t29*t30*t47*t52*(1.0/2.0)-t15*t45*t48*t53*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0)-t30*t34*t45*t48*t52*(1.0/4.0));
  A0[1][1] = -lambda1*(t77+t80+t95+t98+t103-t15*t29*2.0-t20*t24*t30-t15*t24*t43+t20*t30*t50*t52+t15*t43*t50*t55-t29*t30*t52*t53-t15*t45*t53*t58-t29*t34*t41*t96*(3.0/4.0)-t15*t34*t59*t102*(3.0/4.0)+t24*t30*t43*t52*t55*(1.0/2.0)-t30*t34*t45*t52*t58*(1.0/2.0));
  A0[1][2] = -lambda1*(t99+t100+t101+t106+t107+t108+t114+t116-t29*t30*t52*t64*(1.0/2.0)-t29*t30*t53*t63*(1.0/2.0)-t15*t45*t53*t68*(1.0/2.0)-t15*t45*t58*t64*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0)-t30*t34*t45*t52*t68*(1.0/4.0)-t30*t34*t45*t58*t63*(1.0/4.0)-t15*t34*t58*t59*t68*(3.0/4.0));
  A0[1][3] = t130;
  A0[1][4] = lambda1*(t78+t98+t104+t140-t15*t24*t43-t6*t20*t30*t52*(1.0/2.0)-t6*t15*t43*t55*(1.0/2.0));
  A0[1][5] = lambda1*(t99+t100+t101-t7*t20*t30*t52*(1.0/2.0)-t7*t15*t43*t55*(1.0/2.0));
  A0[1][6] = -lambda1*t160;
  A0[1][7] = -lambda1*t175;
  A0[1][8] = -lambda1*(t115+t117+t118-t7*t29*t30*t52*(1.0/2.0)-t7*t15*t45*t58*(1.0/2.0));
  A0[1][9] = lambda1*(t73+t86+t92-t94+t181+t203-t2*t20*t30*t52*(1.0/2.0)-t2*t15*t43*t55*(1.0/2.0)-t29*t30*t40*t53*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0));
  A0[1][10] = lambda1*(t78-t79+t80+t103+t104-t105+t182+t224+t225-t20*t24*t30-t3*t20*t30*t52*(1.0/2.0)-t3*t15*t43*t55*(1.0/2.0)-t29*t30*t52*t53*(1.0/2.0)-t29*t34*t41*t96*(3.0/4.0));
  A0[1][11] = lambda1*(t106+t107+t108+t184+t240-t4*t20*t30*t52*(1.0/2.0)-t4*t15*t43*t55*(1.0/2.0)-t29*t30*t53*t63*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0)-t30*t34*t45*t58*t63*(1.0/4.0));
  A0[2][0] = -lambda1*(t74+t75+t76+t87+t88+t89+t109-t110+t111-t112-t113-t29*t30*t40*t64*(1.0/2.0)-t29*t30*t47*t63*(1.0/2.0)-t15*t45*t48*t64*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0)-t30*t34*t45*t48*t63*(1.0/4.0));
  A0[2][1] = -lambda1*(t99+t100+t101+t106+t107+t108+t114-t115+t116-t117-t118-t29*t30*t52*t64*(1.0/2.0)-t29*t30*t53*t63*(1.0/2.0)-t15*t45*t58*t64*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0)-t30*t34*t45*t58*t63*(1.0/4.0));
  A0[2][2] = -lambda1*(t77+t80+t95+t121+t123-t15*t29*2.0-t20*t24*t30-t15*t24*t43+t20*t30*t61*t63+t15*t43*t61*t66-t29*t30*t63*t64-t15*t45*t64*t68-t29*t34*t41*t119*(3.0/4.0)-t15*t34*t59*t122*(3.0/4.0)+t24*t30*t43*t63*t66*(1.0/2.0)-t30*t34*t45*t63*t68*(1.0/2.0));
  A0[2][3] = t132;
  A0[2][4] = t142;
  A0[2][5] = lambda1*(t78+t121+t124+t149-t15*t24*t43-t7*t20*t30*t63*(1.0/2.0)-t7*t15*t43*t66*(1.0/2.0));
  A0[2][6] = -lambda1*t163;
  A0[2][7] = -lambda1*t178;
  A0[2][8] = -lambda1*t190;
  A0[2][9] = lambda1*(t76+t89+t111-t113+t193+t205-t2*t20*t30*t63*(1.0/2.0)-t2*t15*t43*t66*(1.0/2.0)-t29*t30*t40*t64*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0));
  A0[2][10] = lambda1*(t101+t108+t116-t118+t194+t227-t3*t20*t30*t63*(1.0/2.0)-t3*t15*t43*t66*(1.0/2.0)-t29*t30*t52*t64*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0));
  A0[2][11] = lambda1*(t78-t79+t80+t123+t124-t125+t195+t243+t244-t20*t24*t30-t4*t20*t30*t63*(1.0/2.0)-t4*t15*t43*t66*(1.0/2.0)-t29*t30*t63*t64*(1.0/2.0)-t29*t34*t41*t119*(3.0/4.0));
  A0[3][0] = t128;
  A0[3][1] = t130;
  A0[3][2] = t132;
  A0[3][3] = lambda1*(-t70+t133+t5*t15*t35*t43);
  A0[3][4] = t144;
  A0[3][5] = t151;
  A0[3][9] = -lambda1*(t78+t82-t136-t2*t15*t35*t43*(1.0/2.0));
  A0[3][10] = t229;
  A0[3][11] = t246;
  A0[4][0] = lambda1*(t71+t72+t73-t139-t6*t20*t30*t40*(1.0/2.0));
  A0[4][1] = lambda1*(t78+t98+t104-t133+t140-t6*t20*t30*t52*(1.0/2.0)-t6*t15*t43*t55*(1.0/2.0));
  A0[4][2] = t142;
  A0[4][3] = t144;
  A0[4][4] = lambda1*(-t98+t133+t6*t15*t43*t55);
  A0[4][5] = t153;
  A0[4][9] = t209;
  A0[4][10] = -lambda1*t230;
  A0[4][11] = t248;
  A0[5][0] = lambda1*(t74+t75+t76-t147-t7*t20*t30*t40*(1.0/2.0));
  A0[5][1] = lambda1*(t99+t100+t101-t148-t7*t20*t30*t52*(1.0/2.0));
  A0[5][2] = lambda1*(t78+t121+t124-t133+t149-t7*t20*t30*t63*(1.0/2.0)-t7*t15*t43*t66*(1.0/2.0));
  A0[5][3] = t151;
  A0[5][4] = t153;
  A0[5][5] = lambda1*(-t121+t133+t7*t15*t43*t66);
  A0[5][9] = t212;
  A0[5][10] = t233;
  A0[5][11] = -lambda1*t249;
  A0[6][0] = -lambda1*t157;
  A0[6][1] = -lambda1*t160;
  A0[6][2] = -lambda1*t163;
  A0[6][6] = -lambda1*(t77-t155+t5*t15*t45*t48);
  A0[6][7] = -lambda1*t179;
  A0[6][8] = -lambda1*t191;
  A0[6][9] = t214;
  A0[6][10] = -lambda1*t234;
  A0[6][11] = -lambda1*t250;
  A0[7][0] = -lambda1*(t91+t93+t94-t172-t6*t29*t30*t40*(1.0/2.0));
  A0[7][1] = -lambda1*t175;
  A0[7][2] = -lambda1*t178;
  A0[7][6] = -lambda1*t179;
  A0[7][7] = -lambda1*(t77-t173+t6*t15*t45*t58);
  A0[7][8] = -lambda1*t192;
  A0[7][9] = -lambda1*t216;
  A0[7][10] = t236;
  A0[7][11] = -lambda1*t251;
  A0[8][0] = -lambda1*(t110+t112+t113-t186-t7*t29*t30*t40*(1.0/2.0));
  A0[8][1] = -lambda1*(t115+t117+t118-t187-t7*t29*t30*t52*(1.0/2.0));
  A0[8][2] = -lambda1*t190;
  A0[8][6] = -lambda1*t191;
  A0[8][7] = -lambda1*t192;
  A0[8][8] = -lambda1*(t77-t188+t7*t15*t45*t68);
  A0[8][9] = -lambda1*t218;
  A0[8][10] = -lambda1*t238;
  A0[8][11] = t253;
  A0[9][0] = t201;
  A0[9][1] = lambda1*(t73+t86+t92-t94+t181-t202+t203-t2*t20*t30*t52*(1.0/2.0)-t29*t30*t40*t53*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0));
  A0[9][2] = lambda1*(t76+t89+t111-t113+t193-t204+t205-t2*t20*t30*t63*(1.0/2.0)-t29*t30*t40*t64*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0));
  A0[9][3] = -lambda1*(t78+t82-t136-t206);
  A0[9][4] = t209;
  A0[9][5] = t212;
  A0[9][6] = t214;
  A0[9][7] = -lambda1*t216;
  A0[9][8] = -lambda1*t218;
  A0[9][9] = -lambda1*(t80+t197-t219-t220-t2*t20*t30*t40+t8*t29*t30*t40);
  A0[9][10] = -lambda1*t239;
  A0[9][11] = -lambda1*t254;
  A0[10][0] = lambda1*(t84+t85+t86-t159+t168+t221-t223-t3*t20*t30*t40*(1.0/2.0)-t29*t30*t47*t52*(1.0/2.0)-t29*t34*t40*t41*t52*(3.0/4.0));
  A0[10][1] = lambda1*(t78-t79+t80+t103+t104-t105+t182-t219+t224+t225-t3*t20*t30*t52*(1.0/2.0)-t3*t15*t43*t55*(1.0/2.0)-t29*t30*t52*t53*(1.0/2.0)-t29*t34*t41*t96*(3.0/4.0));
  A0[10][2] = lambda1*(t101+t108+t116-t118+t194-t226+t227-t3*t20*t30*t63*(1.0/2.0)-t29*t30*t52*t64*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0));
  A0[10][3] = t229;
  A0[10][4] = -lambda1*t230;
  A0[10][5] = t233;
  A0[10][6] = -lambda1*t234;
  A0[10][7] = t236;
  A0[10][8] = -lambda1*t238;
  A0[10][9] = -lambda1*t239;
  A0[10][10] = -lambda1*(t80+t103-t219-t3*t20*t30*t52+t9*t29*t30*t52-t29*t34*t41*t96*(3.0/4.0));
  A0[10][11] = -lambda1*t255;
  A0[11][0] = lambda1*(t87+t88+t89-t162+t170+t222-t241-t4*t20*t30*t40*(1.0/2.0)-t29*t30*t47*t63*(1.0/2.0)-t29*t34*t40*t41*t63*(3.0/4.0));
  A0[11][1] = lambda1*(t106+t107+t108-t177+t184+t240-t242-t4*t20*t30*t52*(1.0/2.0)-t29*t30*t53*t63*(1.0/2.0)-t29*t34*t41*t52*t63*(3.0/4.0));
  A0[11][2] = lambda1*(t78-t79+t80+t123+t124-t125+t195-t219+t243+t244-t4*t20*t30*t63*(1.0/2.0)-t4*t15*t43*t66*(1.0/2.0)-t29*t30*t63*t64*(1.0/2.0)-t29*t34*t41*t119*(3.0/4.0));
  A0[11][3] = t246;
  A0[11][4] = t248;
  A0[11][5] = -lambda1*t249;
  A0[11][6] = -lambda1*t250;
  A0[11][7] = -lambda1*t251;
  A0[11][8] = t253;
  A0[11][9] = -lambda1*t254;
  A0[11][10] = -lambda1*t255;
  A0[11][11] = -lambda1*(t80+t123-t219-t4*t20*t30*t63+t10*t29*t30*t63-t29*t34*t41*t119*(3.0/4.0));
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
print matlab_hess_to_c_code(dog_star_var_dict_bnd, hessian_str)
#print double_me(double_me_str)
#print matlab_cross_grad_to_c_code(grad_str)
#print matlab_cross_hessian_to_c_code(hessian_str)
#print matlab_cross_bnd_hessian_to_c_code(hessian_str)
#print matlab_edge_hessian_to_c_code(hessian_str)