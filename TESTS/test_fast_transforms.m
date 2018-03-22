function test_fast_transforms
% NOTE: This code tests the fast fourier2sph/sph2fourier transforms.

clc

% Test 1:
n = 32;
F = sphones(n, n);
err = norm(F - fourier2sph(sph2fourier(F)))/norm(F);
fprintf('Error: %1.1e\n', err)

% Test 2:
n = 32;
f = spherefun.sphharm(0, 0);
ll = linspace(0, 2*pi, 2*n+1);
tt = linspace(0, pi, n+1);
[ll, tt] = meshgrid(ll, tt);
V1 = f(ll, tt);
F = zeros(n+1, 2*n+1);
F(1, 1) = 1;
G = sph2fourier(F);
V2 = fourierCoeffs2vals(G);
err = norm(V1 - V2)/norm(V1);
fprintf('Error: %1.1e\n', err)

% Test 3:
n = 4;
f = spherefun.sphharm(1, 0);
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);
V1 = f(ll, tt);
F = zeros(n+1, 2*n+1);
F(2, 1) = 1;
G = sph2fourier(F);
V2 = fourierCoeffs2vals(G);
err = norm(V1 - V2)/norm(V1);
fprintf('Error: %1.1e\n', err)

% Test 4:
n = 4;
f = spherefun.sphharm(1, -1);
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);
V1 = f(ll, tt);
F = zeros(n+1, 2*n+1);
F(1, 2) = 1;
G = sph2fourier(F);
V2 = fourierCoeffs2vals(G);
err = norm(V1 - V2)/norm(V1);
fprintf('Error: %1.1e\n', err)

% Test 5:
n = 4;
f = spherefun.sphharm(2, -1);
ll = -pi + (0:2*n)*2*pi/(2*n+1);
tt = (0.5:n+1-0.5)*pi/(n+1);
[ll, tt] = meshgrid(ll, tt);
V1 = f(ll, tt);
F = zeros(n+1, 2*n+1);
F(2, 2) = 1;
G = sph2fourier(F);
V2 = fourierCoeffs2vals(G);
err = norm(V1 - V2)/norm(V1);
fprintf('Error: %1.1e\n', err)

end