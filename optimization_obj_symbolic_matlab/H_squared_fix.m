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


s1 = 0.5*(norm(ex_f)+norm(ex_b));
s2 = 0.5*(norm(ey_f)+norm(ey_b));

delta_xf = ex_f/norm(ex_f);
delta_yf = ey_f/norm(ey_f);
delta_xb = ex_b/norm(ex_b);
delta_yb = ey_b/norm(ey_b);

cos_half_b1 = 0.5*norm(delta_xf+delta_xb);
cos_half_b2 = 0.5*norm(delta_yf+delta_yb);

% s1 = avg_len1, s2 = avg_len2, A_star = A = s1*s2
% A(kx+ky) = s1s2(4cos_half_b1/2s1 + 4cos_half_b2/2s2)
%          = 2s1s2(cos_half_b1/s1 + cos_half_b2/s2)
%          = 2(s2cos_half_b1+s1cos_half_b2)
% AH = 0.5(Akx+ky) = s2cos_half_b1+s1cos_half_b2

IntegratedH = simplify(s2*cos_half_b1+s1*cos_half_b2);
IntegratedKx = simplify(s1*cos_half_b1); % not sure what's the best option here

ccode(IntegratedH ,'file','Integrated_H_E');
ccode(gradient(IntegratedH ),'file','Integrated_H_G');

ccode(IntegratedKx,'file','Kx_integrated_E');
ccode(gradient(IntegratedKx),'file','Kx_integrated_G');
H = hessian(Kx_other);
ccode(H,'file','Kx_reg_H');


%H = hessian(absH);
%ccode(H,'file','H_abs_reg_H');

