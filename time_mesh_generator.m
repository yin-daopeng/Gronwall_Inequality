function [tau,t] = time_mesh_generator(Time, N, opt_mesh, r)
if nargin < 3
    deg_graded = "What is the degree of graded mesh? ";
    r = input(deg_graded);
end
%     [T,X]=meshgrid(cumsum([0,2*(Nt+1-(1:Nt))/Nt/(Nt+1)]),(0:dx:1)); %
%     Zhang Ya-nan

switch opt_mesh
    case 1
        t = Time*linspace(0, 1, N+1).^r;
        tau = diff(t);
    case 2
        tau = randi([1,5],[1,1000])/5/N;
        fprintf('rand max tau %d \n', rats(max(tau), 3))
        t = [0, cumsum(tau)];
    case 3
        t1 = linspace(0, .1, 10+1).^r;
        tau1 = .01*rand(1, N-10);
        t = [t1, t1(end)+cumsum(tau1)];
        tau = diff(t);
    case 4
        tau0 = Time/N; % tau min
        t = [];
        t = [t,tau0*rand];
        while t < Time
            r = tau0*rand;
            t = [t, t(end) + r];
        end
        t(end) = Time;
        tau = diff(t);
    case 5
        e = rand(1,N-15); % Li & Liao 
        tau = [(1:15)/N^r, e/sum(e)];
        t = [0,cumsum(tau)];
    case 33
        disp('rho > 1')
        tau = rand(1,N);
        tau = sort(tau);
        % tau = tau/sum(tau);
        t = [0,cumsum(tau)];
    case 44
        disp('rho < 1')
        tau = rand(1,N);
        tau = sort(tau,"descend");
        % tau = tau/sum(tau);
        t = [0,cumsum(tau)];
    case 22
        deg_graded=[1, 9*rand(1,N-1)];
        tau=zeros(1,N);
        tau(1)=.1;
        for i=2:N
            tau(i)=min(tau(i-1)*deg_graded(i), 20);
        end
        t = [0,cumsum(tau)];
    case 55
        disp('uniform mesh')
        tau = .3*rand*ones(1,N);
        t = [0, cumsum(tau)];
    case 6
        tau1 = rand(1,N);
        tau2 = rand(1,N); tau2 = sort(tau2);
        tau3 = rand(1,N); tau3 = sort(tau3,"descend");
        tau4 = rand*ones(1,N);
        tau = [tau2, tau3, tau2, tau3];
        t = [0, cumsum(tau)];
    case 7
        tau1= .4*sin(3*(1:N)/N*pi) + .41 + 0.0*rand(1,N);
        tau = [tau1,tau1];
        t = [0,cumsum(tau)];
    case 8
        tau= .4*sin(16*(1:N)/N*pi) + .41;
        t = [0,cumsum(tau)];
end
end
