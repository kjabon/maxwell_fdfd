%% Defining constants and domain info
u0 = 4e-7*pi;   % H/m
c = 299792458;  % m/s
Zo = u0*c;      % Impedance of free space
n=.8; %resolution (higher = more)
lambda = 1.2179;
k0 = 2*pi/lambda; 
dx = [1/50/n;1/50/n];
NxIn = [500*n; 500*n];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [60*n,60*n;60*n,60*n];
pmlMag = 1;
numDimensions = 2;

%% Make Domain
myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);

%% Defining/calculating positions
grb = 0.1; wb = 0.15; wr = 0.2; R = 1.2; grr = 0.36;
y2 = grr/2+wr/2+2*R+wr+grb+wb/2;  % Bus center distance from oo center.
y1 = grr/2+wr/2+R;  % Ring center distance from oo center.
midpoint = myDomain.GetMidpoint();
srcLocation(1,1) = 1.5;
srcLocation(2,1) = midpoint(2,1)-y2;
xRange = [-1, 10];
yBus1 = [midpoint(2)-y2-wb/2, midpoint(2)-y2+wb/2];
yBus2 = [midpoint(2)+y2-wb/2,midpoint(2)+y2+wb/2];

%% Adding shapes
myDomain = myDomain.AddRectangleInterpolatedEdges('mu', 3, [xRange ; yBus1]);
myDomain = myDomain.AddRectangleInterpolatedEdges('mu', 3, [xRange ; yBus2]);
myDomain = myDomain.AddRing('mu', 3, [midpoint(1),midpoint(2)+y1], R,wr);
myDomain = myDomain.AddRing('mu', 3, [midpoint(1),midpoint(2)-y1], R,wr);
myDomain.PlotGrid('mu',false);

%% Make operators
D = DifferentialOperator(myDomain);
delSquared = D.Laplacian;
mu_zz = myDomain.GetPropertySubgrid('mu_zz',true); %Not quite accurate. Needs fixing
nSquared = (k0*mu_zz).^2;
A = delSquared + nSquared;
fwhmFactor = 2.3548;
fwhm = wb*1.7;
w0 = fwhm*lambda;

%% Make source
Jz = myDomain.Gaussian(srcLocation, 's',2, 1, w0/fwhmFactor);
b = 1j*k0*Zo*Jz;

%% Solve
x = A\b;

%% Plot
x = reshape(x, myDomain.NxS(1), myDomain.NxS(2)); x = x.';
figure; imagesc(real(x));
axis image; set(gca, 'ydir', 'normal');
colormap(redbluehilight);
fmax = max(abs(caxis)); caxis(fmax*[-1 1]); colorbar;