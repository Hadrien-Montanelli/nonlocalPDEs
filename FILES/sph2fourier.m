function G = sph2fourier(F)
%SPH2FOURIER   Convert spherical harmonic coefficients to Fourier coefficients.

loadFastTransforms();

G = jl.call('sph2fourier', F);

end
