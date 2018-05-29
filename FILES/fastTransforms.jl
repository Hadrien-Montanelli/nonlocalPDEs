using FastTransforms

import FastTransforms: leg2cheb, cheb2leg, sph2fourier, fourier2sph

leg2cheb(c::Vector{MxArray}) = [leg2cheb(jvalue(c)) for c in c]'
cheb2leg(c::Vector{MxArray}) = [cheb2leg(jvalue(c)) for c in c]'

sph2fourier(A::Vector{MxArray}) = [sph2fourier(jvalue(A); sketch = :none) for A in A]'
fourier2sph(A::Vector{MxArray}) = [fourier2sph(jvalue(A); sketch = :none) for A in A]'

function fourier_coeffs2vals(A::Matrix)
    P = FastTransforms.plan_synthesis(A)
    B = zero(A)
    A_mul_B!(B, P, A)
end

function fourier_vals2coeffs(A::Matrix)
    P = FastTransforms.plan_analysis(A)
    B = zero(A)
    A_mul_B!(B, P, A)
end
