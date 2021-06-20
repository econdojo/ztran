function   H = CalcH(H_base, H_c, H_s,eta)

global T2 T3 n r

H = H_base + H_c*eta'*T2 + H_s*eta'*T3;

H = H_base + [zeros(n,r) H_c*eta'];

