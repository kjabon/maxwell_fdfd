% Gaussian source

lambda = 1;
k = 2*pi/lambda; 
frequency = 1/lambda;
omega = 2*pi*frequency; %2pi c/lambda == k_0  (careful, c = 1 here)
mu = 1;
eps = 1; %Later these need to be 1d arrays, i.e. an expanded ND matrix

dx = [1/30;1/30];
NxIn = [200; 200];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [35,35;35,35];
pmlMag = 1; 
numDimensions = 2;

myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
%D = DifferentialOperator(myDomain);
%dxSquared = D.dxts*D.dxst;
%dySquared = D.dyts*D.dyst;
%delSquared = dxSquared+dySquared;
%kSquaredMat = diag(sparse(repelem(k*k, size(delSquared,1))));
%A = delSquared + kSquaredMat;

clear xMidpoint;
xMidpoint(1,1) = myDomain.x(round(length(myDomain.x)/2));
xMidpoint(2,1) = myDomain.y(round(length(myDomain.y)/2));

J = myDomain.Gaussian([3;3],2,1,1);
imagesc(J);
axis image;
set(gca, 'ydir', 'normal');