clear all
close all
clc

etaH = 0.15; % NT step size
etaM = 0.05; % NC step size
ttau = 300; % shape parameter
llambda = 2.2; % alpha has weibull distribution with scale parameter lambda and shape parameter kappa
kkapa = 5;
xxi = 10;

K = ones(6) - 0.2 * eye(6);
S = zeros(8,7);
iter = 0;

for i = 1:6
    iter = iter + 1;
    display(['Calculating sensitivity relative to parameter number ' num2str(iter)])
    
    [F, M1, M2, M3, M4, M5, M6] = gen_moments(K(i,1)*etaH, K(i,2)*etaM, K(i,3)*ttau, K(i,4)*llambda, K(i,5)*kkapa, K(i,6)*xxi);
    S(i, :) = [F, M1, M2, M3, M4, M5, M6 + 1836]; 
    % record the objective function and the moments that the model 
    % produces when we increase parameter values by 10%
end
S(7,:) = [0.048, 0.4988, 0.105, 0.6917, 0.0574, 0.608, 1850];
S(8,:) = [0, 0.30, 0.10, 0.60, 0.03, 0.60, 1870];

T = table;
T.var_changed = {'etaH'; 'etaM'; 'tau'; 'lambda'; 'kappa'; 'xi'; 'none'; 'data'};
T.F = S(:,1);
T.new_combo_1880 = S(:,2);
T.new_tech_1880 = S(:,3);
T.new_combo_1930 = S(:,4);
T.new_tech_1930 = S(:,5);
T.peak_reuse = S(:,6);
T.year_peak = S(:,7)



rowLabels = {'$\eta^H$'; '$\eta^M$'; '$\tau$'; '$\lambda$'; '$\kappa$'; '$\xi$'; 'none'; 'data'};
columnLabels = {'Parameter changed'; 'New comb 1880'; 'New tech 1880'; 'New comb 1930'; 'New tech 1930'; 'Peak of reuse'; 'Year of peak'};
matrix2latex(S(:,2:end), 'sensitivity_table_20_new.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels)




