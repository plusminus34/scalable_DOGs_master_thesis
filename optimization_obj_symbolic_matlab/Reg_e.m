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

ex_f = p_xf-p_0;
ey_f = p_yf-p_0;
ex_b = p_xb-p_0;
ey_b = p_yb-p_0;

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

ex_f = p_xf-p_0;
ey_f = p_yf-p_0;
ex_b = p_xb-p_0;
ey_b = p_yb-p_0;

% Difference of squared norm
E_x = simplify((norm(ex_f).^2-norm(ex_b).^2)^2);

ccode(E_x,'file','Reg_E');
ccode(gradient(E_x),'file','Reg_G');
ccode(hessian(E_x),'file','Reg_H');

hess = hessian(E_x);
% seems that if d=abs(norm(ex_f)-norm(ex_b)) than we can use the hessian of
% norm(ex_f).^4+norm(ex_b).^4-(2-d)*norm(ex_f)^2*norm(ex_b)^2) and it is
% tight for this kind of expression!
% since diff will be small this should be reasonable!

diff = sym('diff','real');
E_x2 = simplify(norm(ex_f).^4+norm(ex_b).^4-(2-4*diff)*norm(ex_f)^2*norm(ex_b)^2);

hess_E_x2 = hessian(E_x2);
hess_E_x2 = hess_E_x2(1:9,1:9);
ccode(hess_E_x2,'file','Reg_convex_H');

diff = 0.001;

hess2 = hessian(E_x2);
check = double(subs(hess,[p_0,p_xf,p_xb],[0,0,0 ...
1+diff,0,0,0,0,-1]))
check2 = double(subs(hess2,[p_0,p_xf,p_xb],[0,0,0 ...
1+diff,0,0,0,0,-1]))