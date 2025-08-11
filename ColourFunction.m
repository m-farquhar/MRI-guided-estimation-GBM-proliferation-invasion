function c = ColourFunction(colmat, cmin, cmax, colval)
% Given a colourmap colmat, minimum value in the colormap cmin, maximum
% value in the colormap cmax, determines the color of colval by
% interpolating over RGB space.

[m,~] = size(colmat);

p = polyfit([cmin, cmax], [0,1],1);

idx = polyval(p, colval);
if idx>1
    idx = 1;
elseif idx<0
    idx = 0;
end
x = linspace(0,1,m)';
r = interp1(x(:) ,colmat(:,1), idx);
g = interp1(x(:), colmat(:,2), idx);
b = interp1(x(:), colmat(:,3), idx);

c = [r,g,b];
