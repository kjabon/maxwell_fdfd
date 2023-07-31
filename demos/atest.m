x = 1:50;

for i = 1:length(x)
    bWith1(i) = 1 + 1j*i;
    bNo1(i) = 1j*i;
end
bInvWith1 = 1./bWith1;
bNo1 = 1./bNo1;