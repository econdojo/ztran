function s = SumEta(eta)

global h r
s=0;

for i=0:h-1
    s = s+eta(1+i*r:4+i*r);
end


