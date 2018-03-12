function u = test_nonlocal_ac
% NOTE: This code solves the *nonlocal* Allen-Cahn equation on the sphere with a
% sphharm-based spectral method and ETDRK4. Based on Mikael's fast transform.
% The nonlocal AC equation is: u_t = epsilon^2*L_delta(u) + u - u^3.

% Plotting parameters:
lw = 1; fs = 14; ms = 18;
set(0, 'defaultlinelinewidth', lw)
set(0, 'defaultaxesfontsize', fs)
set(0, 'defaultlinemarkersize', ms)

% Operator parameters:
epsilon = 1e-1;
delta = 1;
alpha = -0.5;

% Discretization parameters:
n = 63;
h = 2e-1;
T = 30;
iterPlot = 5;

% Nonlinearity:
g = @(u) u - u.^3;
c2v = @(u) fourierCoeffs2vals(sph2fourier(u));
v2c = @(u) fourier2sph(fourierVals2coeffs(u));
N = @(u) v2c(g(c2v(u)));
 
% Intial condition:
u0 = randnfunsphere(0.4);

% Get spherical harmonic coefficients of initial condition:
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);
vals = u0(ll, tt);
coeffs = fourierVals2coeffs(vals);
v0 = fourier2sph(coeffs);

% Compute eigenvalues:
lambda = ones(n+1, 1);
for l = 1:n+1
    lambda(l) = evaluateLambda(l-1, delta, alpha);
end
lambda = (1+alpha)*2^(2-alpha)/(delta^2)*lambda;

% Nonlocal Laplace-Beltrami operator: (n+1)(2n+1)x(n+1)(2n+1) diagonal matrix 
% stored as a (n+1)x(2n+1) matrix that acts pointwise.
L = zeros(n+1, 2*n+1);
for i = 1:n+1
    L(i, 1) = lambda(i);
end
for i = 1:n
    for j = 1:n-i+1
        L(i, 2*j) = lambda(i+j);
        L(i, 2*j+1) = lambda(i+j);
    end
end
L = epsilon^2*L;

% Cesaro (C,2) matrix:
A = zeros(n+1, 2*n+1);
for i = 1:n+1
    A(i, 1) = nchoosek(n-i+1+2, n-i+1);
end
for i = 1:n
    for j = 1:n-i+1
        A(i, 2*j) = nchoosek(n-(i+j)+3, n-(i+j)+1);
        A(i, 2*j+1) = A(i, 2*j);
    end
end
A = A/nchoosek(n+2, n);

% ETDRK4 coefficients:
E = exp(h*L); 
E2 = exp(h*L/2);
M = 32; 
r = exp(1i*pi*((1:M) - 0.5)/M);
LR = h*repmat(L(:), 1, M) + repmat(r, (n+1)*(2*n+1), 1);
Q = h*real(mean((exp(LR/2)-1)./LR , 2)); 
f1 = h*real(mean((-4 - LR + exp(LR).*(4 - 3*LR + LR.^2))./LR.^3, 2)); 
f2 = h*real(mean((2 + LR + exp(LR).*(-2 + LR))./LR.^3, 2)); 
f3 = h*real(mean((-4 - 3*LR - LR.^2 + exp(LR).*(4 - LR))./LR.^3, 2)); 
f1 = reshape(f1, n+1, 2*n+1); 
f2 = reshape(f2, n+1, 2*n+1); 
f3 = reshape(f3, n+1, 2*n+1); 
Q = reshape(Q, n+1, 2*n+1);

% Initialize movie:
S = spinopsphere('ac');
pref = spinprefsphere;
pref.iterplot = iterPlot;
pref.grid = 'on';
pref.Clim = [-1 1];
pref.Nplot = 128;
compGrid = getGrid(S, 2*n+2, S.domain);
compGrid = reshapeGrid(S, compGrid);
plotGrid =  getGrid(S, pref.Nplot, S.domain);
plotGrid = reshapeGrid(S, plotGrid);
movieVals = vals2movieVals(vals);
[p, opts] = initializeMovie(S, h, pref, movieVals, compGrid, plotGrid);

% Time-stepping:
itermax = round(T/h);
v = v0;
for iter = 1:itermax
    
    Nv = N(v);
    a = E2.*v + Q.*Nv;
    Na = N(a);
    b = E2.*v + Q.*Na;
    Nb = N(b);
    c = E2.*a + Q.*(2*Nb - Nv);
    Nc = N(c);
    v = E.*v + Nv.*f1 + 2*(Na + Nb).*f2 + Nc.*f3;
    
    % Update movie:
    if ( mod(iter, iterPlot) == 0 )
        
        % Cesaro (C,2) Means:
        vcesaro = A.*v;
        
        % Plot:
        movieVals = vals2movieVals(c2v(vcesaro));
        updateMovie(S, h, p, opts, iter*h, movieVals, compGrid, plotGrid);
    end
    
end

% Output solution from v:
unonloc = sphcoeffs2spherefun(A.*v);

% Create a spinopsphere:
S.lin = @(u) epsilon^2*lap(u);
S.nonlin = @(u) g(u);
S.tspan = [0 T];

% Create initial condition from v0:
S.init = u0;

% Compute and compare with spinsphere result:
pause
uloc = spinsphere(S, 2*n+2, h, pref);
pause
clf
subplot(1, 2, 1), surf(unonloc, 'grid')
caxis([-1 1]), axis equal off
colorbar
title('\textbf{Nonlocal solution}', 'interpreter', 'latex')
subplot(1, 2, 2), surf(uloc, 'grid')
caxis([-1 1]), axis equal off
colorbar
title('\textbf{Local solution}', 'interpreter', 'latex')
err = norm(unonloc - uloc)/norm(uloc);
fprintf('Error: %1.1e\n', err)

% Output both uloc and unonlocal:
u = chebmatrix(unonloc);
u(2, 1) = uloc;

end