% We have 3 edge points (for oscullating plane) and 2 normal points to
% check their angle with their binormals
% Each edge point is defined by two vertices and a fixed parameter t
% The central edge point is p0
% We then need one "inner point" on the fold, v1.
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

% Folded points
fold_v_x = sym('fold_v_x', 'real');
fold_v_y = sym('fold_v_y', 'real');
fold_v_z = sym('fold_v_z', 'real');
fold_v = [fold_v_x,fold_v_y,fold_v_z];

% If there's an isometry energy these should be around the same length
edge = fold_v-ep_0;

% curve binormal vec, not normalized because it's the same from both sides
B = simplify(cross(ep_0-ep_b,ep_f-ep_0));
const = dot(B,edge);

vars = [ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, fold_v_x,fold_v_y,fold_v_z];

ccode(const ,'file','V_fold_vals');
ccode(gradient(const,vars),'file','V_fold_G');

ccode(hessian(const,vars),'file','V_fold_H');