function lambda = evaluateLambda(l, delta, alpha)
%EVALUATELAMBDA   Compute the l-th e-value of the nonlocal Laplacian.

loadFastTransforms();

lambda = jl.call('evaluateLambda', l, delta, alpha);

end
