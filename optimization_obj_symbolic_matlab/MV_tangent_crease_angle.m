ep_b_v1_x = sym('ep_b_v1_x', 'real');
ep_b_v1_y = sym('ep_b_v1_y', 'real');
ep_b_v1_z= sym('ep_b_v1_z', 'real');
ep_b_v1 = [ep_b_v1_x,ep_b_v1_y,ep_b_v1_z];

ep_b_v2_x = sym('ep_b_v2_x', 'real');
ep_b_v2_y = sym('ep_b_v2_y', 'real');
ep_b_v2_z= sym('ep_b_v2_z', 'real');
ep_b_v2 = [ep_b_v2_x,ep_b_v2_y,ep_b_v2_z];

ep_b_t = sym('ep_b_t', 'real'); assume(ep_b_t > 0);

ep_b = ep_b_v1*ep_b_t+(1-ep_b_t)*ep_b_v2;

% Folded points
v1_x = sym('v1_x', 'real');
v1_y = sym('v1_y', 'real');
v1_z = sym('v1_z', 'real');
v1 = [v1_x,v1_y,v1_z];

v2_x = sym('v2_x', 'real');
v2_y = sym('v2_y', 'real');
v2_z = sym('v2_z', 'real');
v2 = [v2_x,v2_y,v2_z];

w1_x = sym('w1_x', 'real');
w1_y = sym('w1_y', 'real');
w1_z = sym('w1_z', 'real');
w1 = [w1_x,w1_y,w1_z];

w2_x = sym('w2_x', 'real');
w2_y = sym('w2_y', 'real');
w2_z = sym('w2_z', 'real');
w2 = [w2_x,w2_y,w2_z];

ep_f_v1_x = sym('ep_f_v1_x', 'real');
ep_f_v1_y = sym('ep_f_v1_y', 'real');
ep_f_v1_z= sym('ep_f_v1_z', 'real');
ep_f_v1 = [ep_f_v1_x,ep_f_v1_y,ep_f_v1_z];

ep_f_v2_x = sym('ep_f_v2_x', 'real');
ep_f_v2_y = sym('ep_f_v2_y', 'real');
ep_f_v2_z= sym('ep_f_v2_z', 'real');
ep_f_v2 = [ep_f_v2_x,ep_f_v2_y,ep_f_v2_z];

ep_f_t = sym('ep_f_t', 'real'); assume(ep_f_t > 0);
ep_f = ep_f_v1*ep_f_t+(1-ep_f_t)*ep_f_v2;

ep_0_t = sym('ep_0_t', 'real'); assume(ep_0_t > 0);
ep_0 = v1*ep_0_t+(1-ep_0_t)*v2;

fold_e_1 = sym('fold_e_1','real');
fold_e_2 = sym('fold_e_2','real');

e1 = (v1-v2)/fold_e_1;
e2 = (w1-w2)/fold_e_2;

l1 = sym('l1','real');
l2 = sym('l2','real');
curve_T = (l2*(ep_0-ep_b)+l1*(ep_f-ep_0))/(l1*l2);

fold_e_crease_angle = sym('fold_e_crease_angle','real');

B = cross(e2,curve_T)/sin(fold_e_crease_angle);

cos_angle = sym('cos_angle','real');

const = dot(e1,B)-cos_angle;

vars = [ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, w1_x, w1_y, w1_z, w2_x, w2_y, w2_z];

ccode(const ,'file','MV_tangent_crease_angle_C');
ccode(gradient(const,vars),'file','MV_tangent_crease_angle_G');
%ccode(hessian(const,vars),'file','FoldTangentsAngleConstraint_H');

%E = const.^2;