%demo2
%2D recovering kx^2 (x direction)

kx = 2*pi/10;%2*pi/1550;


dx = [pi/8, pi/8];
Nx = [50, 40];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [0,0;0,0];
pmlMag = 0; 
numDimensions = 2;
kxvec = kx.*kx*ones(1,Nx(1));

%Create domain and diffOperator instances
myDomain = Domain(dx,Nx,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
%For x in 2D, xvec needs to be repeated for each row...
x = myDomain.GetXSubgrid('s',1);
f = sin(kx*x)';
dxSquared = D.dxts*D.dxst;
g = -(dxSquared*f)./f;


plot(g);
hold on;
plot(kxvec);