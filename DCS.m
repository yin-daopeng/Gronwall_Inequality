
N = 16;
q = zeros(N,1);
for j = 1:N
    q(j) = computing(2^j, -1, 2,.3); % close a constant
end

function q = computing(n,p,q,r)
a = sum((n^r - (1:n-1).^r).^p.*(1:n-1).^q);
if min(p,q) > -1
    b = n^(r*p+q+1);
elseif min(p,q) == -1
    b = n^(r*p+q+1)*(1+log(n));
elseif min(p,q) < -1
    b = n^(max(r*p,(r-1)*p+q));
end
q = a/b;
end