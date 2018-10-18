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

%assume constant length!
cos_a = dot(ex_f,ex_b)/(norm(ex_f)*norm(ex_b));
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

Kx_other = simplify((4/(norm(ex_f)^2+norm(ex_b)^2)).^2*norm(ex_f+ex_b)).^2;
Ky_other = simplify((4/(norm(ey_f)^2+norm(ey_b)^2)).^2*norm(ey_f+ey_b)).^2;
ccode(Kx_other,'file','Kx_reg_E');
ccode(gradient(Kx_other),'file','Kx_reg_G');
H = hessian(Kx_other);
ccode(H,'file','Kx_reg_H');

absH = simplify(0.5*(Kx_other+Ky_other));
ccode(absH ,'file','H_abs_reg_E');
ccode(gradient(absH),'file','H_abs_reg_G');
H = hessian(absH);
ccode(H,'file','H_abs_reg_H');

double(subs(Kx_other,[p_0,p_xf,p_xb],[0,0,0 ...
    1,0,0 ...
    -1,0,0]))

double(subs(gradient(Kx_other),[p_0,p_xf,p_xb],[0,0,0 ...
    1,0,0 ...
    -1,0,0]))

