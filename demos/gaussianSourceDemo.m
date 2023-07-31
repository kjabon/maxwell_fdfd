  % Gaussian sources
tic

lambda = 1;
k = 2*pi/lambda; 
frequency = 1/lambda;
omega = 2*pi*frequency; %2pi c/lambda == k_0  (careful, c = 1 here)
mu = 1;
eps = 1; %Later these need to be 1d arrays, i.e. an expanded ND matrix

dx = [0.005;0.005];
NxIn = [4000; 800];
boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
pmls = [20,20;20,20];
pmlMag = 1; 
numDimensions = 2;

myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
D = DifferentialOperator(myDomain);
dxSquared = D.dxts*D.dxst;
dySquared = D.dyts*D.dyst;
delSquared = dxSquared+dySquared;
kSquaredMat = diag(sparse(repelem(k*k, size(delSquared,1))));
A = delSquared + kSquaredMat;

xMidpoint(1,1) = myDomain.x(round(length(myDomain.x)/2));
xMidpoint(2,1) = myDomain.y(round(length(myDomain.y)/2));

%FWHM = 2.3548 sigma
%Sigma = FWHM/2.3548
%FWHM = f(lambda), e.g. f = 5*lambda, 3*lambda, 1* lambda

fwhmFactor = 2.3548;
fwhm = 2;
w0 = fwhm*lambda;


multiplier = -1i;
J = myDomain.Gaussian(xMidpoint, 's', 1, 1, w0/fwhmFactor);
%J = myDomain.Gaussian(xMidpoint-[0;2.5], 1, 1, w0/fwhmFactor) + multiplier*myDomain.Gaussian(xMidpoint+[0;2.5],1,1,w0/fwhmFactor);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,1,.5);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,-1,.5);
%J = myDomain.Gaussian(xMidpoint-2,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,1i,.5);
%J = myDomain.Gaussian(xMidpoint-1,1,1,.5) + myDomain.Gaussian(xMidpoint+2,1,-1i,.5);

b = (J*omega*mu*1i);

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
    filename = 'demo7.gif';
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

%movie(frames,20);

% subplot(2,1,2);
% imagesc(imag(x)); % imagesc(real(x*exp(1i*pi/2)))
% axis image;
% colorbar;
% set(gca, 'ydir', 'normal');

str = sprintf(': %d, dx: %d, Nx: %d, numPml: %d, pmlMag: %d', lambda, dx(1), NxIn(1), pmls(1), pmlMag);
str = strcat('\lambda' , str);
title(str);
pause(0.1);

% filename = strcat('./output/',str);
% filename = strcat(filename,'.png');
% filename = strrep(filename, ':','');
% filename = strrep(filename, '\','');
% saveas(gcf,filename);

toc
