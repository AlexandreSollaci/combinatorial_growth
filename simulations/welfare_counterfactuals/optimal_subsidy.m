clear
close all
clc

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/welfare_counterfactuals/code/')

etaH = 0.2; % NT step size
etaM = 0.05; % NC step size
etaL = 0; % Refinement step size
ttau = 800; % shape parameter for ideas distribution
llambda = 4; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkappa = 7; % and shape parameter kappa
xxi = 150;   % 1/xi is the fraction of feasible combinations

ggamma = .6; % match to labor share of GDP
epsilon = 2; % from Acemoglu, Akcigit, Bloom and Kerr (2013) - pg 21

rr = 0.05; % interest rate
bbeta = 1/(1+rr); % intertemporal discount factor
LL = 1; % production worker
PPi = ggamma/(1-bbeta); % constant part of the price of innovation

nrofinv = 800; % number of inventors
nroffirms = 1750; % chosen so that avg growth rate is approximately 2%

% define parameters for first 180 years
T_past = 180; % 1836 - 2016
params = v2struct(etaH, etaM, etaL, ttau, llambda, kkappa, xxi, ggamma, epsilon, rr, nrofinv, nroffirms);

% run simulation from 1836 to 2016
result1 = simulate_path(T_past, params);

% sotore simulation results
quality1 = result1.quality;
summat1 = result1.summat;
state1 = result1.state;
Mmat1 = result1.Mmat;
Welfare1 = result1.Welfare;

% define parameters for counterfactuals
T_forward = 50;
subsidy_vals = linspace(.5,1,5);
subs_len = length(subsidy_vals);
subsidy = [0, 0; subsidy_vals', zeros(subs_len,1) ; zeros(subs_len,1), subsidy_vals']; 
    % subsidy for new tech, subsidy for new combo

% initialize matrices to store results
Welfare2 = zeros(1, size(subsidy,1));

% run simulation for counterfactual
for i = 1:size(subsidy,1)
    result2 = simulate_path(T_forward, params, subsidy(i,:), quality1, Mmat1, summat1, state1);
    Welfare2(i) = result2.Welfare; 
end

Welfare = [repmat(Welfare1, 1,  size(subsidy,1)); Welfare2];

Welfare0 = Welfare2(1)*ones(1,subs_len);
WelfareNT = Welfare2(2:subs_len+1);
WelfareNC = Welfare2(subs_len+2:end);

WChangeNT = -(WelfareNT - Welfare0) ./ Welfare0;
WChangeNC = -(WelfareNC - Welfare0) ./ Welfare0;

welfare_gains_pp = figure(1);
plot(subsidy_vals, WChangeNT, '--b', subsidy_vals, WChangeNC, 'r')
legend('New technologies', 'New combination', 'location', 'Northwest')
title('Welfare gains (%) associated with subsidy value')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'welfare_pp.png')

welfare_gains_NT = figure(2);
plot(subsidy_vals, WelfareNT, '--b', subsidy_vals, Welfare0, 'r')
legend('With subsidy', 'Baseline (so subsidy)', 'location', 'Northwest')
title('Welfare gains associated with subsidy for new technologies')
xlabel('Subsidy value')
ylabel('Welfare value')
saveas(gcf, 'welfare_NT.png')

welfare_gains_NC = figure(3);
plot(subsidy_vals, WelfareNC, '--b', subsidy_vals, Welfare0, 'r')
legend('With subsidy', 'Baseline (no subsidy)', 'location', 'Northwest')
title('Welfare gains associated with subsidy for new combinaitons')
xlabel('Subsidy value')
ylabel('Welfare value')
saveas(gcf, 'welfare_NC.png')

%################ Change means of the cost distribution
%
% clear
% close all
% clc

% cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/welfare_counterfactuals/code/')

% etaH = 0.2; % NT step size
% etaM = 0.05; % NC step size
% etaL = 0; % Refinement step size
% ttau = 800; % shape parameter for ideas distribution
% lambda = linspace(1, 5, 5); % alpha (cost) has weibull distribution with scale parameter lambda ...
% kkappa = 7; % and shape parameter kappa
% xxi = 150;   % 1/xi is the fraction of feasible combinations

% ggamma = .6; % match to labor share of GDP
% epsilon = 2; % from Acemoglu, Akcigit, Bloom and Kerr (2013) - pg 21

% rr = 0.05; % interest rate
% bbeta = 1/(1+rr); % intertemporal discount factor
% LL = 1; % production worker
% PPi = ggamma/(1-bbeta); % constant part of the price of innovation

% nrofinv = 800; % number of inventors
% nroffirms = 1750; % chosen so that avg growth rate is approximately 2%

% % define parameters for first 180 years
% T_past = 180; % 1836 - 2016


% % define parameters for counterfactuals
% T_forward = 50;
% subsidy_val = .9;
% lambda_len = length(lambda);
% subsidy = [0, 0; subsidy_val, 0; 0, subsidy_val]; 
%     % subsidy for new tech, subsidy for new combo

% % initialize matrices to store results
% Welfare2 = zeros(lambda_len, size(subsidy,1));

% for j = 1:lambda_len
%     llambda = lambda(j);
% 	params = v2struct(etaH, etaM, etaL, ttau, llambda, kkappa, xxi, ggamma, epsilon, rr, nrofinv, nroffirms);
% 	% run simulation from 1836 to 2016
% 	result1 = simulate_path(T_past, params);

% 	% sotore simulation results
% 	quality1 = result1.quality;
% 	summat1 = result1.summat;
% 	state1 = result1.state;
% 	Mmat1 = result1.Mmat;
% 	Welfare1 = result1.Welfare;

% 	% run simulation for counterfactual
% 	for i = 1:size(subsidy,1)
% 	    result2 = simulate_path(T_forward, params, subsidy(i,:), quality1, Mmat1, summat1, state1);
% 	    Welfare2(j,i) = result2.Welfare; 
% 	end
% end

% WGainNT = -(Welfare2(:,2) - Welfare2(:,1) ) ./ Welfare2(:,1) ;
% WGainNC = -(Welfare2(:,3) - Welfare2(:,1) ) ./ Welfare2(:,1) ;

% figure(1)
% subplot(2,2,1)
% plot(lambda, Welfare2(:,1), 'r', lambda, Welfare2(:,2), 'b', lambda, Welfare2(:,3), 'k');
% legend('no subsidy', 'new tech subsidy', 'new combo subsidy', 'location', 'Southwest')
% title('Welfare of th economy under different cost distributions')
% xlabel('mean of cost distribution'); ylabel('Welfare')
% saveas(gcf, 'welfare_diff_cost.png')

% subplot(2,2,2)
% plot(lambda, Welfare2(:,1), 'r', lambda, Welfare2(:,2), '--b')
% legend('no subsidy', 'subsidy', 'location', 'Southwest')
% title('New technology subsidy')
% xlabel('mean of cost distribution'); ylabel('Welfare')

% subplot(2,2,3)
% plot(lambda, Welfare2(:,1), 'r', lambda, Welfare2(:,3), '--b')
% legend('no subsidy', 'subsidy', 'location', 'Southwest')
% title('New combination subsidy')
% xlabel('mean of cost distribution') ; ylabel('Welfare')

% subplot(2,2,4)
% plot(lambda, WGainNT, '--b', lambda, WGainNC, 'r')
% legend('New tech subsidy', 'new comb subsidy','location', 'Southwest')
% title('Welfare gains from subsidies')
% xlabel('mean of cost distribution') ; ylabel('% Welfare gain')
% saveas(gcf, 'welfare_diff_cost.png')