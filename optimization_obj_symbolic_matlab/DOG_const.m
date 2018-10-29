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


lambda1 = sym('lambda1', 'real');
lambda2 = sym('lambda2', 'real');
lambda3 = sym('lambda3', 'real');

%doesn't assume constant length!
%E_12 = simplify(dot(ex_f,ey_f)*norm(ex_b)-dot(ey_f,ex_b)*norm(ex_f));
%E_23 = simplify(dot(ey_f,ex_b)*norm(ey_b)-dot(ey_b,ex_b)*norm(ey_f));
%E_34 = simplify(dot(ex_b,ey_b)*norm(ex_f)-dot(ex_f,ey_b)*norm(ex_b));

E_12 = simplify(dot(ex_f,ey_f)/(norm(ex_f)*norm(ey_f))-dot(ey_f,ex_b)/(norm(ey_f)*norm(ex_b)));
E_23 = simplify(dot(ey_f,ex_b)/(norm(ey_f)*norm(ex_b))-dot(ey_b,ex_b)/(norm(ey_b)*norm(ex_b)));
E_34 = simplify(dot(ex_b,ey_b)/(norm(ex_b)*norm(ey_b))-dot(ex_f,ey_b)/(norm(ex_f)*norm(ey_b)));

%E_12 = simplify((dot(ex_f,ey_f)*norm(ex_b)).^2-(dot(ey_f,ex_b)*norm(ex_f)).^2);
%E_23 = simplify((dot(ey_f,ex_b)*norm(ey_b)).^2-(dot(ey_b,ex_b)*norm(ey_f)).^2);
%E_34 = simplify((dot(ex_b,ey_b)*norm(ex_f)).^2-(dot(ex_f,ey_b)*norm(ex_b)).^2);

E_all = [E_12;E_23;E_34];

J_3_inner = [diff(E_all,p0_x),diff(E_all,p0_y),diff(E_all,p0_z)...
    diff(E_all,pxb_x),diff(E_all,pxb_y),diff(E_all,pxb_z)...
    diff(E_all,pxf_x),diff(E_all,pxf_y),diff(E_all,pxf_z)...
    diff(E_all,pyb_x),diff(E_all,pyb_y),diff(E_all,pyb_z)...
    diff(E_all,pyf_x),diff(E_all,pyf_y),diff(E_all,pyf_z)];
ccode(J_3_inner ,'file','DOG_Jacobian');
H = lambda1*hessian(E_12,symvar(E_all))+lambda2*hessian(E_23,symvar(E_all))+lambda3*hessian(E_34,symvar(E_all));
ccode(H,'file','DOG_Const_Hessian');
H_bnd = lambda1*hessian(E_12);
ccode(H_bnd,'file','DOG_Const_Hessian_bnd');

Jacobian(E_all)

ccode(E_all,'file','DOG_Constraints_vals');

cos_a = dot(ex_f,ey_f)/(norm(ex_f)*norm(ey_f));
grad_cos_a = gradient(cos_a);

grad_cos_a_sub = subs(grad_cos_a,[p_0,p_xf,p_yf,p_xb,p_yb],[1,2,0 ...
    1,1,0 ...
    3,1,0 ...
    -1,4,0 ...
    5,-1,0]);

J_3_inner_side = [diff(E_all,pxb_x),diff(E_all,pxb_y),diff(E_all,pxb_z)...
    diff(E_all,pyb_x),diff(E_all,pyb_y),diff(E_all,pyb_z)];

J_3_inner_sub = subs(J_3_inner,[p_0,p_xf,p_yf,p_xb,p_yb],[0.1,0.1,0.1 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]);

J_3_inner_side_sub =  subs(J_3_inner_side,[p_0,p_xf,p_yf,p_xb,p_yb],[0,0,0 ...
    1,0,0 ...
    0,1,0.1 ...
    -1,0,0 ...
    0,-1,0]);

E_all_sub = subs(E_all,[p_0,p_xf,p_yf,p_xb,p_yb],[0,0,0 ...
    1,0,0 ...
    0,1,0.1 ...
    -1,0,0 ...
    0,-1,0]);

J_3_inner_1_side = [diff(E_all,pyb_x),diff(E_all,pyb_y),diff(E_all,pyb_z)];
J_3_inner_side_1_sub =  subs(J_3_inner_1_side ,[p_0,p_xf,p_yf,p_xb,p_yb],[0,0,0 ...
    1,0,0 ...
    0,1,0 ...
    -1,0,0 ...
    0,-1,0]);


J_bnd_3_inner = [diff(E_12,p0_x),diff(E_12,p0_y),diff(E_12,p0_z)...
    diff(E_12,pxb_x),diff(E_12,pxb_y),diffnn(E_12,pxb_z)...
    diff(E_12,pxf_x),diff(E_12,pxf_y),diff(E_12,pxf_z)...
    diff(E_12,pyb_x),diff(E_12,pyb_y),diff(E_12,pyb_z)...
    diff(E_12,pyf_x),diff(E_12,pyf_y),diff(E_12,pyf_z)];

ccode(J_bnd_3_inner ,'file','DevOrthJacobianBnd3');
ccode(E_12,'file','DevOrthJacobianValsBnd');

E_corner = simplify(dot(ex_f,ey_f));
J_corner = [diff(E_corner,p0_x),diff(E_corner,p0_y),diff(E_corner,p0_z)...
    diff(E_corner,pxb_x),diff(E_corner,pxb_y),diff(E_corner,pxb_z)...
    diff(E_corner,pxf_x),diff(E_corner,pxf_y),diff(E_corner,pxf_z)...
    diff(E_corner,pyb_x),diff(E_corner,pyb_y),diff(E_corner,pyb_z)...
    diff(E_corner,pyf_x),diff(E_corner,pyf_y),diff(E_corner,pyf_z)];

ccode(J_corner ,'file','DevOrthJacobianBndCorner');
