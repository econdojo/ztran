function [beta, P] = CalcBeta(P,H,M,NQN, old_beta)

global h I_states 

MM = (I_states-old_beta * H)*P;
P = M * MM * M' + NQN;
beta = P * H'*inv(H*P*H');



