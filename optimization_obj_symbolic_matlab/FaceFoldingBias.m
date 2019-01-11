f1_p0_x = sym('f1_p0_x', 'real');
f1_p0_y = sym('f1_p0_y', 'real');
f1_p0_z = sym('f1_p0_z', 'real');

f1_p0 = [f1_p0_x,f1_p0_y,f1_p0_z];

f1_p1_x = sym('f1_p1_x', 'real');
f1_p1_y = sym('f1_p1_y', 'real');
f1_p1_z = sym('f1_p1_z', 'real');

f1_p1 = [f1_p1_x,f1_p1_y,f1_p1_z];

f1_p2_x = sym('f1_p2_x', 'real');
f1_p2_y = sym('f1_p2_y', 'real');
f1_p2_z = sym('f1_p2_z', 'real');

f1_p2 = [f1_p2_x,f1_p2_y,f1_p2_z];

f1_p12_x = sym('f1_p12_x', 'real');
f1_p12_y = sym('f1_p12_y', 'real');
f1_p12_z = sym('f1_p12_z', 'real');

f1_p12 = [f1_p12_x,f1_p12_y,f1_p12_z];

f2_p0_x = sym('f2_p0_x', 'real');
f2_p0_y = sym('f2_p0_y', 'real');
f2_p0_z = sym('f2_p0_z', 'real');

f2_p0 = [f2_p0_x,f2_p0_y,f2_p0_z];

f2_p1_x = sym('f2_p1_x', 'real');
f2_p1_y = sym('f2_p1_y', 'real');
f2_p1_z = sym('f2_p1_z', 'real');

f2_p1 = [f2_p1_x,f2_p1_y,f2_p1_z];

f2_p2_x = sym('f2_p2_x', 'real');
f2_p2_y = sym('f2_p2_y', 'real');
f2_p2_z = sym('f2_p2_z', 'real');

f2_p2 = [f2_p2_x,f2_p2_y,f2_p2_z];

f2_p12_x = sym('f2_p12_x', 'real');
f2_p12_y = sym('f2_p12_y', 'real');
f2_p12_z = sym('f2_p12_z', 'real');

f2_p12 = [f2_p12_x,f2_p12_y,f2_p12_z];


f1_diag1 = f1_p12-f1_p0;
f1_diag2 = f1_p1-f1_p2;

f2_diag1 = f2_p12-f2_p0;
f2_diag2 = f2_p1-f2_p2;

n1 = cross(f1_diag1,f1_diag2);
n2 = cross(f2_diag1,f2_diag2);

n_dot_prod = n1(1)*n2(1)+n1(2)*n2(2)+n1(3)*n2(3);% don't know why but n1.dot(n2) returns an error
obj = simplify(-1*n_dot_prod);
g = gradient(obj);
H = hessian(obj);

fake_obj_n1 = cross(f1_p12,-f1_p2);
fake_obj_n2 = cross(f2_p12,-f2_p2);
fake_obj_dot_prod =  fake_obj_n1(1)*fake_obj_n2(1)+fake_obj_n1(2)*fake_obj_n2(2)+fake_obj_n1(3)*fake_obj_n2(3);
fake_H = hessian(fake_obj_dot_prod);
fake_H_sub = subs(fake_H,[f1_p0,f1_p1,f1_p2,f1_p12,f2_p0,f2_p1,f2_p2,f2_p12],[0,0,0,1,0,0,0,1,0,1/sqrt(2),1/sqrt(2),0, 0,0,0,1,0,0,0,1,0,1/sqrt(2),1/sqrt(2),0]);
double(eig(fake_H_sub))