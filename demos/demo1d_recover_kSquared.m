%1D recovering kx^2
kx = 2*pi/10;%2*pi/1550;


dx = 1;
Nx = 5;
boundaryConditions = [BCs.antiSymS,BCs.antiSymT];
pmls = [0,0];
pmlMag = 0;
numDimensions = 1;
kxvec = kx.*kx*ones(1,Nx);

%Create domain and diffOperator instances
myDomain = Domain(dx,Nx,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
x = myDomain.GetXSubgrid('s',1);
f = sin(kx*x-pi/8)';
dxSquared = D.dxts*D.dxst;
g = -(dxSquared*f)./f;
plot(g);
hold on;
plot(kxvec);