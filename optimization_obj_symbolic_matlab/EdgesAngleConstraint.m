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

e1 = v1-v2;
e2 = w1-w2;

cos_angle = sym('cos_angle','real');

const = dot(e1,e2)-cos_angle;

vars = [v1_x, v1_y, v1_z, v2_x, v2_y, v2_z,w1_x, w1_y, w1_z, w2_x, w2_y, w2_z];

ccode(const ,'file','EdgesAngleConstraint_C');
ccode(gradient(const,vars),'file','EdgesAngleConstraint_G');
%ccode(hessian(const,vars),'file','FoldTangentsAngleConstraint_H');

%E = const.^2;