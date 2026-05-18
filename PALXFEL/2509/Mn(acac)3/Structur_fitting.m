
clc; clear; close all;
%% Data
load("Mn_U_prime.mat");
aff_Mn = [11.2519 5.34818 7.36935 0.34373 3.04107 17.4089 2.27703 84.2139 1.05195];
aff_O = [2.95648 13.8964 2.45240 5.91765 1.50510 0.34537 0.78135 34.0811 0.30413];
aff_C = [1.93019 12.7188 1.87812 28.6498 1.57415 0.59645 0.37108 65.0337 0.24637];

heating = load("MeOH_comps.dat");
q_heat = heating(:,1);
heat = q_heat.*heating(:,2:4);
clear heating
heat = interp1(q_heat,heat,q,'linear','extrap');
heat = smoothdata(heat,1,"movmean",10);

f_Mn = aff(aff_Mn,q);
f_O = aff(aff_O,q);
f_C = aff(aff_C,q);

%%
chi2_func = @(r) 10*(U_prime(:,3) + r(13)*heat(:,1) + r(14)*heat(:,2) + r(15)*heat(:,3) - r(16)*debye(r,q,f_Mn,f_O,f_C)  );

lb = [0;0;0;...
      0;0;0;...
      
      1;1;1;...
      1;1;1;...

      0;0;0;...

      0.00001];

ub = [4;4;4;...
      4;4;4;...

      4;4;4;...
      4;4;4;...

      2;2;2;...

      0.001];

x0 = [2;2;2;...
      2;2;2;...

      1;1;1;...
      1;1;1;...
      
      1;1;1;...
      
      0.0005];
%% 
option = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 2e5, ...
    'StepTolerance',1e-6);

problem = createOptimProblem('lsqnonlin', ...
                             'objective', chi2_func, ...
                             'options', option, ...
                             'x0', x0, ...
                             'lb', lb, ...
                             'ub', ub);


ms = MultiStart('Display', 'iter', 'UseParallel',true);
[p_final, fval, exitflag, output, solutions] = run(ms, problem, 1000);


%%

DS = p_final(16)*debye(p_final(1:12),q,f_Mn,f_O,f_C) ; %+ sum(p_final(9:11)'.*heat,2)
figure;
hold on;
plot(q,U_prime(:,3) + sum(p_final(13:15)'.*heat,2),'o');
plot(q,DS);
legend('data','fit');

r = zeros(12,1);
for i = 1 : 6
    d2 = p_final(6+i);
    d1 = p_final(i);
    r(2*i-1) = norm([d1;d2]-[1.4;0]);
    r(2*i) = norm([d1;d2]-[-1.4;0]);
    
end





%%
% function DS = debye(r,q,f1,f2)
% DS = f1.*f2.*(math_sinc(q*r(4)) + math_sinc(q*r(5)) + math_sinc(q*r(3)) ...
%             - 2*math_sinc(q*r(1)) - math_sinc(q*r(2)));
% DS = q.*DS;
% end

function DS = debye(r,q,f1,f2,f3)
DS = acac(r(4),r(10),q,f1,f2,f3) + acac(r(5),r(11),q,f1,f2,f3) + acac(r(12),r(8),q,f1,f2,f3) ...
    -acac(r(1),r(7),q,f1,f2,f3) -acac(r(2),r(8),q,f1,f2,f3) - acac(r(3),r(9),q,f1,f2,f3);
DS = q.*DS;
end

function S = acac(d1,d2,q,fMn,fO,fC)
Mn = [d1;d2];
O1 = [1.3955;0];
O2 = -O1;

C1 = [1.2238; -1.2452];
C2 = [-1.2238; -1.2452];

C3 = [0; -1.8936];

C4 = [3.25; -1.60];
C5 = [-3.25; -1.60];


ring1 = 2*fMn.*fO.*(math_sinc(q*norm(Mn-O1)) + math_sinc(q*norm(Mn-O2)));
ring2 = 2*fMn.*fC.*(math_sinc(q*norm(Mn-C1)) + math_sinc(q*norm(Mn-C2)));
ring3 = 2*fMn.*fC.*(math_sinc(q*norm(Mn-C3)) + math_sinc(q*norm(Mn-C4)) + math_sinc(q*norm(Mn-C5)));

S = ring1 + ring2 + ring3;

end