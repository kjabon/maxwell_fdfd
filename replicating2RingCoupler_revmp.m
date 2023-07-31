% Ken's replica of the 2-ring EC568 example

close all; clearvars;
%addpath(genpath('..'));

% Defining constants and domain info
u0 = 4e-7*pi; c = 299792458; Zo = u0*c;  % H/m, m/s, ohms (free space impedance)

n=.5;               % Resolution (higher = more)
dx = [1/50/n;1/50/n];
NxIn = [500*n; 500*n];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [60*n,60*n;60*n,60*n];  pmlMag = 1;
numDimensions = 2;

% Make Domain
myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
%myPML = PML(width, maxconductivity, order);
%myDomain = Domain([xmin:dx/2:xmax],[PEC PEC],[myPML myPML]); % should not have to
%specify S or T -- grid as given will determine that, need to just specify
%starting pixel as S or T I guess..  How do we best tell it a reference
%point so it knows which pixels are S and T? (maybe AA, AB, BA and BB grids..)

% Defining/calculating positions
grb = 0.1; wb = 0.15; wr = 0.2; R = 1.2; grr = 0.36;
y2 = grr/2+wr/2+2*R+wr+grb+wb/2;    % Bus center distance from oo center.
y1 = grr/2+wr/2+R;                  % Ring center distance from oo center.
midpoint = myDomain.GetMidpoint();
srcLocation(1,1) = 1.5;
srcLocation(2,1) = midpoint(2)-y2;
xRange = [-1, 10];
yBus1 = [midpoint(2)-y2-wb/2, midpoint(2)-y2+wb/2];
yBus2 = [midpoint(2)+y2-wb/2,midpoint(2)+y2+wb/2];

% Adding shapes
myDomain = myDomain.AddRectangleInterpolatedEdges('eps', 3^2, [xRange ; yBus1]);
myDomain = myDomain.AddRectangleInterpolatedEdges('eps', 3^2, [xRange ; yBus2]);
myDomain = myDomain.AddRing('eps', 3^2, [midpoint(1),midpoint(2)+y1], R,wr);
myDomain = myDomain.AddRing('eps', 3^2, [midpoint(1),midpoint(2)-y1], R,wr);
myDomain.PlotGrid('eps',false);

% Make operators
D = DifferentialOperator(myDomain);
delSquared = D.Laplacian;

lamvec = linspace(1.216,1.220,10);


for lambda = lamvec(2)
    %lambda = 1.2179;
    k0 = 2*pi/lambda; 
    disp(sprintf('Wavelength %g um...', lambda));

    eps_zz = myDomain.GetPropertySubgrid('eps_zz',true); %Not quite accurate. Needs fixing
    k0sqnSquared = k0^2*eps_zz;
    A = delSquared + k0sqnSquared;

    % Make source
    fwhmFactor = 2.3548; fwhm = wb*1.7; w0 = fwhm*lambda;
    Jz = myDomain.Gaussian(srcLocation, 's', 2, 1, w0/fwhmFactor);
    b = 1j*k0*Zo*Jz;

    % Solve
    x = A\b;

    % Plot
    x = reshape(x, myDomain.NxS(1), myDomain.NxS(2)); x = x.';
    figure; imagesc(real(x));
    axis image; set(gca, 'ydir', 'normal');
    colormap(redbluehilight);
    fmax = max(abs(caxis)); caxis(fmax*[-1 1]); colorbar;
    drawnow;
end

Ezmax = fmax; %max(abs(Ezmat(:)));
figure;
for mm = 1:24
    wt = (mm-1)/24 * 2*pi;
%    imagesc(real(xEz), real(yEz), real(Ezmat.'*exp(j*wt)));
    imagesc(real(x * exp(j*wt)));
    set(gca,'ydir','normal'); axis image; colorbar; colormap(redbluehilight);
    caxis([-1 1]*Ezmax);
    xlabel('Position x (\mum)'); ylabel('Position y (\mum)'); title('E_z radiation given J_z source (2D)');
    set(gcf,'Color',[1 1 1]);
    
    M(mm) = getframe;
end
movie(M,5)
