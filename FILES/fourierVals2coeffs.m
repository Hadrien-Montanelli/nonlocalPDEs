function G = fourierVals2coeffs(F)
%FOURIERCOEFFS2VALS   Convert values on a grid of equispaced angles on the 
%sphere to Fourier coefficients.

loadFastTransforms();

G = jl.call('fourier_vals2coeffs', F);

end