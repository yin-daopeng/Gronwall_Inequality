clear
clc

Alpha = .6;
Beta = .3;
Nt = 2.^(6:9);
N = length(Nt);
dim =1;

% N*(meshtype +1)*(normtype x2)
Etable = zeros(N,7,2);
Etable(:,1,:)=repmat(Nt',[1,1,2]); % error table

%================== setting degree of grading =============================
r0=(2-Alpha)/(1+Beta-Alpha);
R = [1, r0, r0 + 1];

for r = 1:length(R)
    for n = 1:N
        [~,er] = solver_Dalhquist(Nt(n), Alpha, Beta, R(r));
        Etable(n,2*r,1) = er(2);
        Etable(n,2*r,2) = er(end);
    end
end

for i = 1:2
    for r = 1:length(R)
        for n = 2:N
            Etable(n,2*r+1,i) = - log(Etable(n,2*r,i)/Etable(n-1,2*r,i))...
                                 /log(Etable(n,  1,i)/Etable(n-1,  1,i));
        end
        if i ==1
        fprintf( '[iCo: %2.3f/%2.3f]', R(r)*Beta,min([2-Alpha,R(r)*(1+Beta-Alpha)]));
        end
    end
end
fprintf('\n')

fprintf('\n ================== Error ============================== \n');
for n=1:N
    fprintf('%5d  &  %1.3e & %1.3f &  %1.3e & %1.3f &  %1.3e & %1.3f \\\\ \n',...
        Etable(n,:,1));
end
fprintf('\n --------------------------------------------------------- \n');
for n=1:N
    fprintf('%5d  &  %1.3e & %1.3f &  %1.3e & %1.3f &  %1.3e & %1.3f \\\\ \n',...
        Etable(n,:,2));
end


