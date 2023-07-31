function J = redbluehilight(m)
%REDBLUEHILIGHT    Color scheme with +ve values blue, -ve red, zero white.
%   REDBLUEHILIGHT(M), a variant of JET(M), is an M-by-3 matrix containing
%   the default colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red.
%   JET, with no arguments, is the same length as the current colormap.
%   Use COLORMAP(JET).
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   March 18, 2009, Milos Popovic.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
n = floor(m/4);
if (mod(m,2) == 1) ODD = 1; else ODD = 0; end;
uL = [0:1:n-1]/(n-1 +ODD);
uR = fliplr(uL);
if (ODD) uC = 1; else uC = []; end;
uo = ones(1,n);

% r = [(uo+uL)/2 uo uC uR uo*0].';
% g = [uo*0 uL uC uR uo*0].';
% b = [uo*0 uL uC uo (uo+uR)/2].';
r = [uo uo uC uR uo*0].';
g = [uR uL uC uR uL].';
b = [uo*0 uL uC uo uo].';
J = [r g b];

% To test the colormap try:
% imagesc(sin([0:0.02:100].'/4)*(cos([0:0.02:100]/4).*exp(-(([0:0.02:100]/4)-10).^2/50) ))
% colormap(redbluehilight)

% [X,Y]=meshgrid(-5:.02:5,-5:.02:5);
% Z = exp(-(sqrt(X.^2+Y.^2)-3).^2/0.25) .* exp(i*10*atan2(Y,X));
% imagesc(real(Z)); caxis([-1 1]); axis image; colorbar; colormap(redbluehilight);
% modeanim(X(1,:).',Y(:,1),Z, 2, 1);
