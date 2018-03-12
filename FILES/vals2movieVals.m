function vvv = vals2movieVals(vv)
%VALS2MOVIEVALS   Interpolate values used by Mikael's fast transform to a grid 
%that can be used by spinopsphere movie functions.

n = size(vv, 1) - 1;

% Grid used by Mikael's fast transform:
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);

% Grid used by spherefun:
lll = trigpts(2*n + 2, [-pi pi])';
ttt = linspace(0, pi, n+1);
[lll, ttt] = meshgrid(lll, ttt);

% Interpolate:
vvv = interp2(ll, tt, vv, lll, ttt, 'spline');

% Double-up grid to pass to spinopsphere movie functions:
vvv = [flipud(vvv); vvv];

end