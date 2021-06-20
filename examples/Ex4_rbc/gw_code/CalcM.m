function   M = CalcM(L,eta,beta,H,M,MW);

global T1
global T2
global T3
global T7
global I_states

M1 = T3*((L*eta'+(I_states-beta * H)*M*T1)*T2+beta*H*M*T7);
    
M =   [ MW;
        M1];
