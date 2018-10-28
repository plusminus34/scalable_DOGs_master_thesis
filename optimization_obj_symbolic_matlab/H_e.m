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

delta_1 = ex_f/norm(ex_f);
delta_2 = ey_f/norm(ey_f);
delta_1b = ex_b/norm(ex_b);
delta_2b = ey_b/norm(ey_b);

%assume constant length!
cos_x = dot(delta_1,delta_1b);
cos_y = dot(delta_2,delta_2b);

alpha = acos(cos_x);
beta = acos(cos_y);

% new H
half_k_x_new = simplify(norm(delta_1+delta_1b)/(norm(ex_f)+norm(ex_b)));
half_k_y_new = simplify(norm(delta_2+delta_2b)/(norm(ey_f)+norm(ey_b)));

% other formulation
Half_K_x2 = simplify(2*(cos(alpha/2)/(norm(ex_f)+norm(ex_b))));
Half_K_y2 = simplify(2*(cos(beta/2)/(norm(ey_f)+norm(ey_b))));

H_old = simplify(half_k_x+half_k_y);
H_new_1 = simplify(half_k_x_new+half_k_y_new);
H_new_2 = simplify(Half_K_x2+Half_K_y2);

E_all = simplify(H_new_1^2);
E_bnd = simplify(half_k_x_new^2);
%E_code = ccode(E_all);

ccode(E_all,'file','H_E');
ccode(E_bnd,'file','H_E_bnd');
ccode(gradient(E_all),'file','H_G');
ccode(gradient(half_k_x_new),'file','H_K_1_G');
ccode(gradient(half_k_x_new+half_k_y_new),'file','H_K_1_plus_K2_G');
ccode(gradient(E_bnd),'file','H_G_bnd');



ccode(hessian(E_all),'file','H_Hess');
ccode(hessian(E_bnd),'file','H_Hess');

Kx_other = (4/(norm(ex_f)^2+norm(ex_b)^2))*norm(ex_f+ex_b);
Ky_other = (4/(norm(ey_f)^2+norm(ey_b)^2))*norm(ey_f+ey_b);

E_new = 0.5*( Kx_other.^2+Ky_other.^2);
E_new = simplify(E_new);
E_other = (4/(norm(ex_f)^2+norm(ex_b)^2)*()

H_old_val = double(subs(H_old,[p_0,p_xf,p_yf,p_xb,p_yb],[0.1,0.1,0.1 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]));

H_new_val = double(subs(H_new_1,[p_0,p_xf,p_yf,p_xb,p_yb],[0.1,0.1,0.1 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]));

H_2_new_val= double(subs(H_new_2,[p_0,p_xf,p_yf,p_xb,p_yb],[0.1,0.1,0.1 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]));

bla = gradient(E_bnd);
g_check = double(subs(bla,[p_0,p_xf,p_yf,p_xb,p_yb],[0,0,0 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]));



