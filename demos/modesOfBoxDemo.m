%demo3
%2D modes of cavity, with neumann boundaries on left and right, and zero
%boundaries on top and bottom

close all;
mu = 1;
eps = 1; %Later these need to be 1d arrays, i.e. an expanded ND matrix
dx = [pi/32, pi/32];
NxIn = [200, 300];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [0,0;0,0];
pmlMag = 0; 
numDimensions = 2; 

%Create domain and diffOperator instances
myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
dxSquared = D.dxts*D.dxst;
dySquared = D.dyts*D.dyst;
delSquared = dxSquared+dySquared;
A = 1*(1/mu)*(1/eps)*delSquared;

[V,eigVals] = eigs(A,6,'smallestabs');
%plot(V(:,1));
Nx = myDomain.NxS(1);
Ny = myDomain.NxS(2);
V1 = reshape(V(:,1),Nx(1), Ny);
V1 = V1.';
V2 = reshape(V(:,2),Nx(1), Ny);
V2 = V2.';
V3 = reshape(V(:,3),Nx(1), Ny);
V3 = V3.';
V4 = reshape(V(:,4),Nx(1), Ny);
V4 = V4.';
V5 = reshape(V(:,5),Nx(1), Ny);
V5 = V5.';
V6 = reshape(V(:,6),Nx(1), Ny);
V6 = V6.';

imagesc(V1)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');
figure
imagesc(V2)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');
figure
imagesc(V3)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');
figure
imagesc(V4)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');
figure
imagesc(V5)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');
figure
imagesc(V6)
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');

