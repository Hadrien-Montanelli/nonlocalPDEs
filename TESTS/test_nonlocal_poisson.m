function u = test_nonlocal_poisson
% NOTE: This code solves the *nonlocal* Poisson equation on the sphere with a
% sphharm-based spectral method. Comparison with local solution using DFS.

clc, clf, drawnow

% Operator parameters:
delta = 1;
alpha = -0.5;

% Discretization parameters:
n = 130;

% Right-hand side:
ff = @(x,y,z) -(exp(-30*((y+sqrt(3)/2).^2+x.^2+(z-1/2).^2)) + exp(-50*z.^2));
f = spherefun(ff);

% Get spherical harmonic coefficients of right-hand side:
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);
vals = f(ll, tt);
coeffs = fourierVals2coeffs(vals);
F = fourier2sph(coeffs);

% Compute eigenvalues:
lambda = ones(n+1, 1);
for l = 1:n+1
    lambda(l) = evaluateLambda(l-1, delta, alpha);
end
lambda = (1+alpha)*2^(2-alpha)/(delta^2)*lambda;

% Inverse of nonlocal Laplace-Beltrami operator: (n+1)(2n+1)x(n+1)(2n+1) 
% diagonal matrix stored as a (n+1)x(2n+1) matrix that acts pointwise.
Linv = zeros(n+1, 2*n+1);
for i = 1:n+1
    Linv(i, 1) = 1/lambda(i);
end
for i = 1:n
    for j = 1:n-i+1
        Linv(i, 2*j) = 1/lambda(i+j);
        Linv(i, 2*j+1) = 1/lambda(i+j);
    end
end
Linv(1, 1) = 1;

% Solve:
U = Linv.*F;

% Output solution from U:
unonloc = sphcoeffs2spherefun(U);

% Compute local result:
uloc = spherefun.poisson(f, mean2(f), 2*n+2);

% Plot:
subplot(1, 3, 1), surf(f), axis off
title('\textbf{Right-hand side}', 'interpreter', 'latex')
subplot(1, 3, 2), surf(unonloc), caxis([-0.18 -0.05]), axis off
title('\textbf{Nonlocal solution}', 'interpreter', 'latex')
subplot(1, 3, 3), surf(uloc), caxis([-0.18 -0.05]), axis off
title('\textbf{Local solution}', 'interpreter', 'latex')
colormap('jet')

% Some quantities:
fprintf('Mean of f: %1.4f\n', mean2(f))
fprintf('Mean of nonlocal u: %1.4f\n', mean2(uloc))
fprintf('Mean of local u: %1.4f\n', mean2(unonloc))
err = norm(unonloc - uloc)/norm(uloc);
fprintf('Error: %1.1e\n', err)

% Output both uloc and unonlocal:
u = chebmatrix(unonloc);
u(2, 1) = uloc;

end