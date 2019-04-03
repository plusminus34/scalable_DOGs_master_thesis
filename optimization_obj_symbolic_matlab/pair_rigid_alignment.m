p1_x = sym('p1_x', 'real');
p1_y = sym('p1_y', 'real');
p1_z = sym('p1_z', 'real');

p1 = [p1_x,p1_y,p1_z];

p2_x = sym('p2_x', 'real');
p2_y = sym('p2_y', 'real');
p2_z = sym('p2_z', 'real');

p2 = [p2_x,p2_y,p2_z];

syms r11 r12 r13 r21 r22 r23 r31 r32 r33;
R = [r11,r12,r13;r21,r22,r23;r31,r32,r33];
syms t1 t2 t3
T = [t1,t2,t3];

% Equal up to rigid motion
const = T+p1*R-p2;
E = sum(const.^2);

vars = [p1_x, p1_y, p1_z, p2_x, p2_y, p2_z];

ccode(E,'file','pair_rigid_alignment_E');
ccode(gradient(E,vars),'file','pair_rigid_alignment_G');
ccode(hessian(E,vars),'file','pair_rigid_alignment_H');