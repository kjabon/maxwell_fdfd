%demo3
%2D recovering ky^2 (y direction)

ky = 2*pi/10;%2*pi/1550;


dx = [pi/16, pi/16];
Nx = [50, 40];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [0,0;0,0];
pmlMag = 0; 
numDimensions = 2;
kyvec = ky.*ky*ones(1,Nx(2)*Nx(1));

%Create domain and diffOperator instances
myDomain = Domain(dx,Nx,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
%For x in 2D, xvec needs to be repeated for each row...
y = myDomain.GetXSubgrid('s',2);
f = sin(ky*y)';
dySquared = D.dyts*D.dyst;
g = -(dySquared*f)./f;
plot(g);
hold on;
plot(kyvec);