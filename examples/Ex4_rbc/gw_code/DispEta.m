function DispEta(n)

% display n orders of the hierachy at period t
global r
global h
global eta

if n > h-1
    n=h-1;
end

ETA = eta(1:r);
for i = 1:n
    ETA = [ETA eta(i*r+1:(i+1)*r)];
end

ETA
