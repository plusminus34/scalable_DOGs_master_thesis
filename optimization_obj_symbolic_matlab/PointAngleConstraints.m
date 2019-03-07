v1_x = sym('v1_x', 'real');
v1_y = sym('v1_y', 'real');
v1_z = sym('v1_z', 'real');
v1 = [v1_x,v1_y,v1_z];

v2_x = sym('v2_x', 'real');
v2_y = sym('v2_y', 'real');
v2_z = sym('v2_z', 'real');
v2 = [v2_x,v2_y,v2_z];

center_x = sym('center_x', 'real');
center_y  = sym('center_y', 'real');
center_z  = sym('center_z', 'real');
center = [center_x,center_y,center_z];

e1 = v1-center;
e2 = center-v2;

cos_angle = sym('cos_angle','real');

const = dot(e1,e2)-cos_angle;

vars = [v1_x, v1_y, v1_z, v2_x, v2_y, v2_z,w1_x, w1_y, w1_z, w2_x, w2_y, w2_z];

ccode(const ,'file','PointAngleConstraints_C');
ccode(gradient(const,vars),'file','PointAngleConstraints_G');
%ccode(hessian(const,vars),'file','FoldTangentsAngleConstraint_H');

%E = const.^2;