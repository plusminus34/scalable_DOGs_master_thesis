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

pxff_x = sym('pxff_x', 'real');
pxff_y = sym('pxff_y', 'real');
pxff_z = sym('pxff_z', 'real');

p_xff = [pxff_x,pxff_y,pxff_z];

ex_f = p_xf-p_0;
ex_b = p_xb-p_0;
ex_ff = p_xff-p_xf;


len_ex_b = sym('len_ex_b');
len_ex_f = sym('len_ex_f');
len_ex_ff = sym('len_ex_ff');

assume(len_ex_b > 0)
assume(len_ex_f > 0)
assume(len_ex_ff > 0)

% Kx*N integrated is linear if we know the lengths
KxN1 = 2*(ex_f/len_ex_f+ex_b/len_ex_b);
KxN2 = 2*(ex_ff/len_ex_ff+ex_f/len_ex_f);
diffKN = KxN1-KxN2;
E = simplify(norm(diffKN).^2);
ccode(E ,'file','fairing_E');
vars = [p_0,p_xb,p_xf,p_xff];
ccode(gradient(E,vars),'file','fairing_G');
ccode(hessian(E, vars),'file','fairing_H');