p0_x = sym('p0_x', 'real');
p0_y = sym('p0_y', 'real');
p0_z = sym('p0_z', 'real');

p_0 = [p0_x,p0_y,p0_z];

pxf_x = sym('pxf_x', 'real');
pxf_y = sym('pxf_y', 'real');
pxf_z = sym('pxf_z', 'real');

p_xf = [pxf_x,pxf_y,pxf_z];

pxb_x = sym('pxb_x', 'real');
pxb_y = sym('pxb_y', 'real');
pxb_z = sym('pxb_z', 'real');

p_xb = [pxb_x,pxb_y,pxb_z];

pyf_x = sym('pyf_x', 'real');
pyf_y = sym('pyf_y', 'real');
pyf_z = sym('pyf_z', 'real');

p_yf = [pyf_x,pyf_y,pyf_z];

pyb_x = sym('pyb_x', 'real');
pyb_y = sym('pyb_y', 'real');
pyb_z = sym('pyb_z', 'real');

p_yb = [pyb_x,pyb_y,pyb_z];

% distance from squared length (l^0 should be the squared length!)
diag1 = simplify(norm(p_yf-p_xf)^2);
diag2 = simplify(norm(p_xb-p_yf)^2);
diag3 = simplify(norm(p_yb-p_xb)^2);
diag4 = simplify(norm(p_xf-p_yb)^2);

% If it would be constraints, we would need only3 but since its an energy
% there's no reason not to make it symmetric
const_1 = diag1-diag2;
const_2 = diag2-diag3;
const_3 = diag3-diag4;
const_4 = diag4-diag1;

E_all = simplify(const_1^2+const_2^2+const_3^2+const_4^2);

ccode(E_all,'file','EqualDiag_E');
ccode(simplify(gradient(E_all)),'file','EqualDiag_G');
ccode(hessian(E_all),'file','EqualDiag_H');

E_bnd = simplify(const_1^2);
ccode(E_bnd,'file','EqualDiagBnd3_E');
ccode(simplify(gradient(E_bnd)),'file','EqualDiagBnd3_G');