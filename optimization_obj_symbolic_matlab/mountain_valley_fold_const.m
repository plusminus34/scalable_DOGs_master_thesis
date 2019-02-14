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

% Fold edge point
ep_0_v1_x = sym('ep_0_v1_x', 'real');
ep_0_v1_y = sym('ep_0_v1_y', 'real');
ep_0_v1_z = sym('ep_0_v1_z', 'real');
ep_0_v1 = [ep_0_v1_x,ep_0_v1_y,ep_0_v1_z];

ep_0_v2_x = sym('ep_0_v2_x', 'real');
ep_0_v2_y = sym('ep_0_v2_y', 'real');
ep_0_v2_z = sym('ep_0_v2_z', 'real');
ep_0_v2 = [ep_0_v2_x,ep_0_v2_y,ep_0_v2_z];

ep_0_t = sym('ep_0_t', 'real'); assume(ep_0_t > 0);
ep_0 = ep_0_t*ep_0_v1+(1-ep_0_t)*ep_0_v2;

% Folded point
folded_v_x = sym('folded_v_x', 'real');
folded_v_y = sym('folded_v_y', 'real');
folded_v_z = sym('folded_v_z', 'real');
folded_v = [folded_v_x,folded_v_y,folded_v_z];

lambda = sym('lambda', 'real');

% If there's an isometry energy these should be around the same length
edge = folded_v-ep_0;

% curve binormal vec, not normalized because it's the same from both sides
B = simplify(cross(ep_0-ep_b,ep_f-ep_0));
const = dot(B,edge);

vars = [ep_0_v1_x,ep_0_v1_y,ep_0_v1_z,ep_0_v2_x,ep_0_v2_y,ep_0_v2_z,ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z,folded_v_x,folded_v_y,folded_v_z];

ccode(const ,'file','V_fold_vals');
ccode(gradient(const,vars),'file','V_fold_G');

ccode(lambda*hessian(const,vars),'file','V_fold_H');

ccode(lambda*hessian(log(const),vars),'file','V_fold_E_H');