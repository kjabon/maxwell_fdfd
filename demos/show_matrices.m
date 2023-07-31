clc
disp('antiSymS, s->t, and antiSymT,, t->s')
disp(full(DantiSymS.dxst) )   %0 both ends as special values
disp('symS, s->t, and symT, t->s')
disp(full(DsymS.dxst) )       %sqrt(2) both ends as special values (- on - diagonal)


%full(DantiSymT.dxts)
%full(DsymT.dxts)


disp('antiSymS, t->s, and antiSymT, s->t')
disp(full(DantiSymS.dxts))
disp('symS, t->s, and symT, s->t')
disp(full(DsymS.dxts))

%full(DantiSymT.dxst)
%full(DsymT.dxst)

for