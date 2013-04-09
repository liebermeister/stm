function [C,Ceq] = elf2(x,N,epsilon)

eval(default('epsilon','10^-5'));
[nm,nr] = size(N);
C       = (diag(x(1:nr)) * N' * x(nr+1:end))+epsilon;
Ceq     = [];
