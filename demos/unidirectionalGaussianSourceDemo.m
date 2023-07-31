  % Gaussian sources
  addpath(genpath('..'));
ns = [.2 .3 .4 .5 .6 .7 .8 .9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0];% 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3];
%ns = [1 2 3 4 5 6 7 8 9 10];
ns=1;
nIndices = 1:length(ns);
clear backMax;
clear launchMax;
clear backRatio;
figure;
for nIndex = nIndices
    
lambda = 1;
k = 2*pi/lambda; 
frequency = 1/lambda;
omega = 2*pi*frequency; %2pi c/lambda == k_0  (careful, c = 1 here)
mu = 1;
eps = 1; %Later these need to be 1d arrays, i.e. an expanded ND matrix
n = ns(nIndex);
disp(n);

dx = [1/30/n;1/30/n];
NxIn = [200*n; 200*n];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [50*n,50*n;50*n,50*n];
pmlMag = 1;
numDimensions = 2;

myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
dxSquared = D.dxts*D.dxst;
dySquared = D.dyts*D.dyst;
delSquared = dxSquared+dySquared;

jOmegaEps = diag(sparse(repelem(1i*omega*eps, size(delSquared,1))));

A = -(omega*mu*1i)^-1*delSquared + jOmegaEps;

xMidpoint(1,1) = myDomain.x(round(length(myDomain.x)/2));
xMidpoint(2,1) = myDomain.y(round(length(myDomain.y)/2));
syvec = GetXSubvector(myDomain, 's', 2);
xMidpoint(2,1) = syvec(length(syvec)/2);

%FWHM = 2.3548 sigma
%Sigma = FWHM/2.3548
%FWHM = f(lambda), e.g. f = 5*lambda, 3*lambda, 1* lambda

fwhmFactor = 2.3548;
fwhm = 2;
w0 = fwhm*lambda;

multiplier = -1i;
Jz = myDomain.Gaussian(xMidpoint, 's',1, 1, w0/fwhmFactor);

%calculate the phase shift in terms of one dx/2 step
phaseStep = dx(2)/2/lambda*2*pi;
phaseStepMultiplier = exp(1i*phaseStep);


%Mx = myDomain.Gaussian(xMidpoint - [0;dx(2)/2], 'u', 1, 1, w0/fwhmFactor);
Mx = 1/2*myDomain.Gaussian(xMidpoint-[0;dx(2)/2], 'u', 1, 1, w0/fwhmFactor)  +...
    1/2*myDomain.Gaussian(xMidpoint + [0;dx(2)/2], 'u', 1, 1, w0/fwhmFactor);
%J = myDomain.Gaussian(xMidpoint-[0;2.5], 1, 1, w0/fwhmFactor) + multiplier*myDomain.Gaussian(xMidpoint+[0;2.5],1,1,w0/fwhmFactor);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,1,.5);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,-1,.5);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,1i,.5);
%J = myDomain.Gaussian(xMidpoint-1,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,-1i,.5);
%3.3000    3.3167    3.3333



b = Jz+D.dyts*(omega*mu*1i)^-1*Mx;

x = A\b;
x = reshape(x, myDomain.NxS(1), myDomain.NxS(2));
x = x.';

imagesc(real(x));
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');

makeGif = false;
if makeGif
    N = 50;
    fmax = max(abs(x(:)));
    filename = 'demo8.gif';
    h = figure;
    for k = 1:N %To stop 1 frame short
        wt = (k-1)/N * 2*pi; % do one full cycle
        fieldwt = x * exp(-1i*wt);
        imagesc(real(fieldwt));
        axis image; 
        set(gca, 'ydir', 'normal');
        colormap(redbluehilight);
        caxis(fmax*[-1 1]);

        drawnow;
        frame = getframe(h); 
        image = frame2im(frame);
        [indexedImg, colorMap] = rgb2ind(image, 256);

        if k == 1
            imwrite(indexedImg, colorMap, filename, 'gif', 'Loopcount', inf,'DelayTime',0);
        else
            imwrite(indexedImg, colorMap, filename, 'gif', 'WriteMode', 'append','DelayTime',0);
        end
    end
end

str = sprintf(': %d, dx: %d, Nx: %d, numPml: %d, pmlMag: %d', lambda, dx(1), NxIn(1), pmls(1), pmlMag);
str = strcat('\lambda' , str);
title(str);
pause(0.1);

backMax(nIndex) = max(abs(real(x(1:round(NxIn(1)/2-5),round(NxIn(2)/2)))));
launchMax(nIndex) = max(abs(real(x(1:round(NxIn(1)/2+5),round(NxIn(2)/2)))));

end
% backRatio = backMax./launchMax;
% figure;
% plot(ns,backRatio);
% xlabel('n (dx = dx_0/n)');
% ylabel('Maximum backward/Maximum forward amplitude');

