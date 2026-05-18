
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
scale_factor = [1; 1; 1; 1; 1e-1; 1; 1e-1; 1; 1e-1; 1e-1; 1e-1; 1e-4];
chi2_func = @(r_scaled) 10*(U_prime(:,3) ...
    + (r_scaled(9)*scale_factor(9)) * heat(:,1) ...
    + (r_scaled(10)*scale_factor(10)) * heat(:,2) ...
    + (r_scaled(11)*scale_factor(11)) * heat(:,3) ...
    - (r_scaled(12)*scale_factor(12)) * debye(r_scaled .* scale_factor, q, f_Mn, f_O, f_C) );

lb = [0.5;0.5;...
      0.5;0.5;...
      0;0;...
      0;0;...
      0;0;0;...
      0.00001];

ub = [3;3;...
      3;3;...
      pi/2;0;...
      pi/2;0;...
      2;2;2;...
      0.001];

x0 = [2;2;...
      2;2;...
      0;0;...
      0;0;...
      1;1;1;...
      0.0005];

lb_scaled = lb ./ scale_factor;
ub_scaled = ub ./ scale_factor;
x0_scaled = x0 ./ scale_factor;
%% 
option = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'FiniteDifferenceType', 'central', ...
    'StepTolerance', 1e-8, ...
    'FunctionTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 5e5, ...
    'MaxIterations', 2000);

problem = createOptimProblem('lsqnonlin', ...
                             'objective', chi2_func, ...
                             'options', option, ...
                             'x0', x0_scaled, ...
                             'lb', lb_scaled, ...
                             'ub', ub_scaled);


ms = MultiStart('Display', 'iter', 'UseParallel',true);
[p_final, fval, exitflag, output, solutions] = run(ms, problem, 1000);

p_final = p_final .* scale_factor;
%%

DS = p_final(12)*debye(p_final(1:8),q,f_Mn,f_O,f_C) ; %+ sum(p_final(9:11)'.*heat,2)
figure;
hold on;
plot(q,U_prime(:,3) + sum(p_final(9:11)'.*heat,2),'o');
plot(q,DS);
legend('data','fit');

r = zeros(12,1);
for i = 1 : 2
    theta = p_final(3+2*i);
    d = [p_final(2*i-1);p_final(2*i)];
    r(6*i-5) = norm([cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;d(1)]-[1.4;0]);
    r(6*i-4) = norm([cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;d(1)]-[-1.4;0]);
    r(6*i-3) = norm([cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;d(1)]-[1.4;0]);
    r(6*i-2) = norm([cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;d(1)]-[-1.4;0]);
    r(6*i-1) = norm([cos(0), -sin(0); sin(0), cos(0)]*[0;d(2)]-[-1.4;0]);
    r(6*i) = norm([cos(0), -sin(0); sin(0), cos(0)]*[0;d(2)]-[-1.4;0]);

end




%%
% function DS = debye(r,q,f1,f2)
% DS = f1.*f2.*(math_sinc(q*r(4)) + math_sinc(q*r(5)) + math_sinc(q*r(3)) ...
%             - 2*math_sinc(q*r(1)) - math_sinc(q*r(2)));
% DS = q.*DS;
% end

function DS = debye(r,q,f1,f2,f3)
DS = 2*acac(r(3),r(7),q,f1,f2,f3) + acac(r(4),r(8),q,f1,f2,f3) ...
    -2*acac(r(1),r(5),q,f1,f2,f3) - acac(r(2),r(6),q,f1,f2,f3);
DS = q.*DS;
end

function S = acac(d,theta,q,fMn,fO,fC)
Mn = [cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;d];
O1 = [1.4; 0];
O2 = -O1;

C1 = [1.23; -1.27];
C2 = [-1.23; 1.27];

C3 = [0; -1.91];

C4 = [2.50; -2.05];
C5 = [-2.5; -2.05];


ring1 = 2*fMn.*fO.*(math_sinc(q*norm(Mn-O1)) + math_sinc(q*norm(Mn-O2)));
ring2 = 2*fMn.*fC.*(math_sinc(q*norm(Mn-C1)) + math_sinc(q*norm(Mn-C2)));
ring3 = 2*fMn.*fC.*(math_sinc(q*norm(Mn-C3)) + math_sinc(q*norm(Mn-C4)) + math_sinc(q*norm(Mn-C5)));

S = ring1 + ring2 + ring3;

end