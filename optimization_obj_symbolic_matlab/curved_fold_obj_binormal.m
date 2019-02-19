% We have 3 edge points (for oscullating plane) and 2 normal points to
% check their angle with their binormals
% Each edge point is defined by two vertices and a fixed parameter t


% Edge points (also need the first one since there is a fold between
% v1,v2..)

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

lambda = sym('lambda','real');

% If there's an isometry energy these should be around the same length
e1 = v1-v2;
e2 = w1-w2;

% curve binormal vec, not normalized because it's the same from both sides
B = simplify(cross(ep_0-ep_b,ep_f-ep_0));

ccode(B,'file','curved_fold_obj_binormal_B_fixed');

% We want <e1,B> = - <e2,B> or <e1,B>+<e2,B> = 0
curve_fold_const = simplify(dot(e1,B) + dot(e2,B));

ccode(curve_fold_const ,'file','FoldingBinormalBiasConstraint_const');
ccode(gradient(curve_fold_const,vars) ,'file','FoldingBinormalBiasConstraint_gradient');
ccode(hessian(lambda*curve_fold_const,vars) ,'file','FoldingBinormalBiasConstraint_lambdaHessian');

E = simplify(curve_fold_const.^2);
% list of variables and their order (the 't' variables are seen as
% fixed parameters and we should  not calculate theirgradient/hessian
% entrires)
vars = [ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z,w1_x, w1_y, w1_z, w2_x, w2_y, w2_z];

ccode(E ,'file','curved_fold_obj_binormal_E');
ccode(gradient(E,vars),'file','curved_fold_obj_binormal_G');

% get a linearized dot(e1,B)+dot(e2,B)

% Not clear how to simplify the hessian for now... so maybe now stay
% without an hessian
H = hessian(E,vars);
ccode(H,'file','curved_fold_obj_binormal_H');

%subH = subs(H,[ep_0_v1_x, ep_0_v1_y, ep_0_v1_z, ep_0_v2_x,ep_0_v2_y,ep_0_v2_z, ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x,ep_f_v2_y, ep_f_v2_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, ep_0_t,ep_b_t,ep_f_t]  ...
%    ,[0,0,0, 0,1,0,-1,0,0,-1,1,0, 1,0,0,1,1,0  0,0,0, 0,1,0,  0.5,0.4,0.5 ]);


x0 = sym('x0',size(vars));
%B_fixed_x = sym('B_fixed_x', 'real');
%B_fixed_y = sym('B_fixed_y', 'real');
%B_fixed_z= sym('B_fixed_z', 'real');
%B_fixed = [B_fixed_x,B_fixed_y,B_fixed_z];
B_fixed = taylor(B,vars,'ExpansionPoint',x0,'Order',1);
E_simplified = simplify((dot(e1+e2,B_fixed)).^2);
H_simp = hessian(E_simplified,vars);
H_simp = subs(H_simp, x0,vars);
ccode(H_simp,'file','curved_fold_obj_binormal_H_simp');
%subH_simp = subs(H_simp,[ ep_0_v1_x, ep_0_v1_y, ep_0_v1_z, ep_0_v2_x, ep_0_v2_y, ep_0_v2_z, ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x,ep_f_v2_y, ep_f_v2_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, ep_0_t,ep_b_t,ep_f_t] ,[0,0,0,0,1,0, -1,0,0,-1,1,0, 1,0,0,1,1,0  0,0,0, 0,1,0.1,  0.5,0.4,0.5 ]);

% Another way to simplify: Find a linear approximation of B
% B has only quadratic coefficients, but we can fix only v1,v2 and at that
% case B is linear in the other coefficient. We can then represent B = B0 +
% linear coefficients
% Which will give a PSD quadratic objective with some hessian that has more
% information than just the 2 vertices.. (for instance states how the other
% two vertices should move as well)
% Last way: Fix ep_0 instead of v1,v2. In that case we again let the other
% parts move (the hessian should be the same?). 
% In any case place B = B_fixed + linear_approx


% Another idea, fix the tangent, and just use the principle normal (as a
% fixed linear combination of the edges). Maybe this makes sense as the
% act of folding the tangent itself is the rotation axis?
% Other idea, somehow fix the angle between the tangent and the principle
% normal or tangent and edges
taylor(B, vars_without_v1_v2, 'Order', 2)