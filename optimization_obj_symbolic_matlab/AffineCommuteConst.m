syms r1_11 r1_12 r1_13 r1_21 r1_22 r1_23 r1_31 r1_32 r1_33;
R1 = [r1_11,r1_12,r1_13;r1_21,r1_22,r1_23;r1_31,r1_32,r1_33];
syms t1_1 t1_2 t1_3
T1 = [t1_1,t1_2,t1_3];

syms r2_11 r2_12 r2_13 r2_21 r2_22 r2_23 r2_31 r2_32 r2_33;
R2 = [r2_11,r2_12,r2_13;r2_21,r2_22,r2_23;r2_31,r2_32,r2_33];
syms t2_1 t2_2 t2_3
T2 = [t2_1,t2_2,t2_3];

assume(R1,'real')
assume(T1,'real')
assume(R2,'real')
assume(T2,'real')

const_R = R1*R2-R2*R1; const_R = const_R(:);
const_T = T2*R1+T1-T1*R2-T2;
const = [const_R;const_T'];
vars = [t1_1,t1_2,t1_3,r1_11,r1_12,r1_13,r1_21,r1_22,r1_23,r1_31,r1_32,r1_33,t2_1,t2_2,t2_3,r2_11,r2_12,r2_13,r2_21,r2_22,r2_23,r2_31,r2_32,r2_33];

ccode(const ,'file','affine_commute_C');
ccode(jacobian(const,vars),'file','affine_commute_G');