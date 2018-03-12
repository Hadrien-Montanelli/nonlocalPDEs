function u = sphcoeffs2spherefun(v)
%SPHCOEFFS2SPHEREFUN   Make a SPHEREFUN object from spherical harmonic coeffs.

n = size(v, 1) - 1;

% Get Fourier coefficients from spherical coeffs:
coeffs = sph2fourier(v);

% Get values from Fourier coefficients:
vv = fourierCoeffs2vals(coeffs);

% Grid used by Mikael's fast transform:
ll = -pi + (0:2*n)*2*pi/(2*n+1); % 2n+1 points
tt = (0.5:n+1-0.5)*pi/(n+1); % n+1 points
[ll, tt] = meshgrid(ll, tt); 

% Grid used by spherefun:
lll = trigpts(2*n + 2, [-pi pi])'; % need even number so I chose 2n+2
ttt = linspace(0, pi, n + 1); % n+1 points is OK
[lll, ttt] = meshgrid(lll, ttt); 

% Interpolate on spherefun sgrid:
vvv = interp2(ll, tt, vv, lll, ttt, 'spline'); 

% Create spherefun from values:
u = spherefun(vvv);

end