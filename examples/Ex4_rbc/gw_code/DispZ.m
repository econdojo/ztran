function DispZ(n,t)

% display n orders of the hierachy at period t
global r
global h
global Z_s ;

if n > h
    n=h;
end

Output = Z_s(1:r,t);
for i = 1:n
    Output = [Output Z_s(i*r+1:(i+1)*r,t)];
end

Output
