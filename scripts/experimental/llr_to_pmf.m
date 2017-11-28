N = 6;

L = sort(abs(randn(1,N/2)));

L = [-fliplr(L), L]
exp(L)

A = -eye(N) + fliplr(diag(exp(L)))

A = [A; ones(1,size(A,2))]
rank(A)
inv_A = pinv(A)
p = ( A\[zeros(N,1); 1]).'

L - log(p./fliplr(p))