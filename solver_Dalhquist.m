
function [uh, er] = solver_Dalhquist(Nt, Alpha, Beta, r)
T = 1;
[tau,t] = time_mesh_generator(T, T*Nt, 1, r);

Kappa = 1; % \kappa 
const = 0;


%--------------------------------------------------------------------------
uh = zeros(1,length(t));
u = @(t) t.^Beta/gamma(1+Beta)+const; % exact solution
f = @(t) t.^(Beta-Alpha)/gamma(Beta-Alpha+1) - Kappa*(t.^Beta/gamma(1+Beta)+const);
%----------------------------------------
uu   = u(t);
ff   = f(t);
uh(1) = uu(1);

%======================= discrete convolution kernels =====================
[k,n]=meshgrid(1:length(tau));
if Alpha == 1
    dck = diag(1./tau);
else
    dck =  ((t(1+n) - t(k)  ).^(1-Alpha) ...
        -   (t(1+n) - t(k+1)).^(1-Alpha))./(t(k+1) - t(k))/gamma(2 - Alpha);
    dck = tril(dck); % asymptotically compatible
end

his = 0;

for n = 1: length(tau)    
    M = dck(n, n)  - Kappa;
    if n > 1; his = (diff(uh(1:n),1,2)*dck(n, 1:n-1)'); end % the history parts;
    v = dck(n, n)*uh(n) - his + ff(n+1);
    uh(n+1) = M\v;
end
er = abs(uu-uh);


