function N = CalcN(beta,H,N_old)

global r
global T3
global T8

N = [   eye(r)
        T3*beta*H*N_old*T8];
