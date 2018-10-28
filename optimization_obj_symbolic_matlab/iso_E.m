p0_x = sym('p0_x', 'real');
p0_y = sym('p0_y', 'real');
p0_z = sym('p0_z', 'real');

p_0 = [p0_x,p0_y,p0_z];

pxf_x = sym('pxf_x', 'real');
pxf_y = sym('pxf_y', 'real');
pxf_z = sym('pxf_z', 'real');

p_xf = [pxf_x,pxf_y,pxf_z];

l0 = sym('l0', 'real');

% distance from squared length (l^0 should be the squared length!)
l = simplify(norm(p_xf-p_0).^2);
assume(l0 > 0)
assume(l > 0)
E_all = simplify((l-l0)^2);

ccode(E_all,'file','Iso_E');
ccode(gradient(E_all),'file','Iso_G');
ccode(hessian(E_all),'file','Iso_H');