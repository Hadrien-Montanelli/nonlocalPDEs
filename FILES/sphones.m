function A = sphones(m, n)
%SPHONES   Create a matrix of spherical harnomic coefficients all equal to 1.

A = zeros(m, 2*n-1);
for i = 1:m
    A(i,1) = 1.0;
end
for j = 1:n-1
    for i = 1:m-j
        A(i,2*j) = 1.0;
        A(i,2*j+1) = 1.0;
    end
end

end
