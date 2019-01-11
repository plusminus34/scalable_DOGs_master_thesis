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


E_12 = simplify(dot(ex_f,ey_f)-dot(ey_f,ex_b));
E_23 = simplify(dot(ey_f,ex_b)-dot(ey_b,ex_b));
E_34 = simplify(dot(ex_b,ey_b)-dot(ex_f,ey_b));

E_all = [E_12,E_23,E_34];
E = simplify(norm(E_all).^2);
E_bnd4 = simplify(E_12.^2+E_23.^2);
E_bnd3 = simplify(E_12.^2);

ccode(E ,'file','DOG_quad_E');
ccode(simplify(gradient(E)),'file','DOG_quad_G');

ccode(E_bnd4,'file','DOG_quad_bnd4_E');
ccode(simplify(gradient(E_bnd4)),'file','DOG_quad_bnd4_G');

ccode(E_bnd3,'file','DOG_quad_bnd3_E');
ccode(simplify(gradient(E_bnd3)),'file','DOG_quad_bnd3_G');