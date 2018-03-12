function F = fourier2sph(G)
%FOURIER2SPH   Convert Fourier coefficients to spherical harmonic coefficients.

loadFastTransforms();

F = jl.call('fourier2sph', G);

end
