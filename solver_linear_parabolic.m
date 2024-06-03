function [uh,H1e,L2e,Linfe] = solver_linear_parabolic(Nt, Nx, Alpha, Beta, r)
% c_{\alpha} \partial_{t}^{\alpha}u - \epsilon \Delta u = \kappa u + f(x,t)
T = 1;
[tau,t] = time_mesh_generator(T, T*Nt, 1, r);

x = linspace(-pi, pi, Nx+1);

Calpha = 1; 
Eps = 1; % \epsilon 
Kap = 1; % \kappa 
const  = 1;


%--------------------------------------------------------------------------
uh = zeros(Nx+1,length(t));
uex = @(x,t)   sin(x).*(t.^Beta/gamma(1+Beta)+const); % exact solution
phi = @(x,t)   cos(x).*(t.^Beta/gamma(1+Beta)+const); % Neumann bdr;
src = @(x,t)   sin(x).*(Calpha*t.^(Beta-Alpha)/gamma(Beta-Alpha+1) ...
                      - Kap*(t.^Beta/gamma(1+Beta)+const) ...
                      + Eps*(t.^Beta/gamma(1+Beta)+const));

%----------------------------------------
[tt,xx] = meshgrid(t,x);
u   = uex(xx, tt);
F   = src(xx, tt);
Phi = phi(xx, tt);
%---------------------- initial & boundary condition ----------------------
uh(:,1) = u(:,1);
uh([1,end],:) = u([1,end],:);
%==========================================================================
h = max(diff(x));
D = toeplitz([-2,1,zeros(1,Nx-3)]); % discrete Laplacian 

chose_bdr = 1;
switch chose_bdr
    case 1 % 'Dirichlet'
        Lh = D/h^2;
        F([2,end-1], :) = F([2,end-1], :) + Eps/h^2*u([1,end],:);
    case 2 % 'Neumann'
        Lh(  1,    1:2  ) = -[ 1, -1]/h^2;
        Lh(end,end-1:end) = -[-1,  1]/h^2; % [-3, 4 , -1]/dx/2;
        F([1,end], :)= 0*F([1,end], :) - 0*[2/h;2/h] + [-1; 1].*Phi([1,end],:)/h;
%     case 3 % 'Robin'
end

%======================= discrete convolution kernels =====================
[k,n]=meshgrid(1:length(tau));
if Alpha == 1
    dck = diag(1./tau);
else
    dck =  ((t(1+n) - t(k)  ).^(1-Alpha) ...
        -   (t(1+n) - t(k+1)).^(1-Alpha))./(t(k+1) - t(k))/gamma(2 - Alpha);
    dck = tril(dck); % asymptotically compatible
end

E = eye(size(Lh));
his = zeros(Nx-1,1);
uh_in = uh(2:end-1,:);

for n = 1: length(tau)    
    M = Calpha*dck(n, n)*E - Eps*Lh - Kap*E;
    if n > 1; his = Calpha*(diff(uh_in(:, 1:n),1,2)*dck(n, 1:n-1)'); end % the history parts;
    v = Calpha*dck(n, n)*uh_in(:, n) - his + F(2:end-1,n+1);
    uh_in(:, n+1) = M\v;
end

uh(2:end-1,2:end) = uh_in(:,2:end);

%% post-processing
L2e  = sqrt(h*sum((u-uh).^2));
H1e =  sqrt(L2e.^2 + sum((diff(u-uh)/h).^2*h));
Linfe = max(abs(u-uh));

figure(1)
subplot(1,3,1)
plot(t,L2e, 'b.-',LineWidth=.5); 
title('error distribution')
xlabel('time')
ylabel('$L^2$ error', Interpreter='latex')

subplot(1,3,2)
plot(t,H1e, 'b.-',LineWidth=.5); 
title('error distribution')
xlabel('time')
ylabel('$H^1$ error',Interpreter='latex')

subplot(1,3,3)
plot(t,Linfe, 'b.-',LineWidth=.5); 
title('error distribution')
xlabel('time')
ylabel('$L^\infty$ error', Interpreter='latex')



