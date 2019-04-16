% We have 3 edge points (for oscullating plane), ep_0,ep_b,ep_f
% We also have two additional "inner folded points" for folding. One of
% them might be the vertex 'v' but it doesn't matter which
% Each edge point is defined by two vertices and a fixed parameter t
% The central edge point is at t*ep_0_v1 + (1-t)ep_0_v2
% The vertex 'v' is the one we fold
% This function is the dot product of v1-p0 with the binormal B
% So the condition that its positive means its a valley fold
% To have a mountain constraint, flip the sign of the function
% (value, gradient and hessian)

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

% Folded points
v1_x = sym('v1_x', 'real');
v1_y = sym('v1_y', 'real');
v1_z = sym('v1_z', 'real');
v1 = [v1_x,v1_y,v1_z];

v2_x = sym('v2_x', 'real');
v2_y = sym('v2_y', 'real');
v2_z = sym('v2_z', 'real');
v2 = [v2_x,v2_y,v2_z];

ep_0_t = sym('ep_0_t', 'real'); assume(ep_0_t > 0);
ep_0 = ep_0_t*v1+(1-ep_0_t)*v2;

e1 = v1-v2;
B = simplify(cross(ep_0-ep_b,ep_f-ep_0));

alpha = sym('alpha', 'real'); assume(ep_0_t > 0);

% We want <e1,B> = - <e2,B> or <e1,B>+<e2,B> = 0
%curve_fold_const = simplify(tanh(alpha*dot(e1,B)) + tanh(alpha*dot(e2,B)));
mv_const = 0.5 + 0.5*tanh(alpha*dot(e1,B));

vars = [ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z];

ccode(mv_const  ,'file','MV_fold_const');
ccode(gradient(mv_const,vars),'file','MV_fold_G');