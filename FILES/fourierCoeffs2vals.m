function G = fourierCoeffs2vals(F)
%FOURIERCOEFFS2VALS   Convert Fourier coefficients to values on a grid of
%equispaced angles on the sphere.

loadFastTransforms();

G = jl.call('fourier_coeffs2vals', F);

end