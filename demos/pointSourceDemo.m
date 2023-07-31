%function  demo6(multiplier1,multiplier2,filename)%(lambda,dx,NxIn,pmls,pmlMag)

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
D = DifferentialOperator(myDomain);
dxSquared = D.dxts*D.dxst;
dySquared = D.dyts*D.dyst;
delSquared = dxSquared+dySquared;
kSquaredMat = diag(sparse(repelem(k*k, size(delSquared,1))));
A = delSquared + kSquaredMat;

xMidpoint(1,1) = myDomain.x(round(length(myDomain.x)/2));
xMidpoint(2,1) = myDomain.y(round(length(myDomain.y)/2));

multiplier1 = 1;

J = multiplier1*myDomain.Delta(xMidpoint,'s');%+multiplier2*myDomain.Delta(xMidpoint+[2;0]);
b = (J*omega*mu*1i);
%b = exp(1i).';

x = A\b;
x = reshape(x, myDomain.NxS(1), myDomain.NxS(2));
x = x.';

% subplot(2,1,1)
imagesc(real(x));
axis image;
colormap(redbluehilight);
fmax = max(abs(caxis));
caxis(fmax*[-1 1]); % Detect the min and max field value, and make the colorbar symmetric about the zero field so that it is white, and +ve and -ve are blue and red respectively.
colorbar;
set(gca, 'ydir', 'normal');

makeGIF = false;
if makeGIF
N = 50;
fmax = max(abs(x(:)));
filename = 'demo6.gif';
filename = strcat(filename,'.gif');
filename = char(filename);
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

% filename = strcat('./output/',str);
% filename = strcat(filename,'.png');
% filename = strrep(filename, ':','');
% filename = strrep(filename, '\','');
% saveas(gcf,filename);

%end
