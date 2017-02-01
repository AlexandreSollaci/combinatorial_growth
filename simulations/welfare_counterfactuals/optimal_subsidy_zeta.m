clear
close all
clc

%cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/welfare_counterfactuals/')

etaH = 0.3; % NT step size
etaM = 0.07; % NC step size
etaL = 0; % Refinement step size
ttau = 500; % shape parameter for ideas distribution
llambda = 1; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkappa = 5; % and shape parameter kappa
xxi = 200;   % 1/xi is the fraction of feasible combinations

zeta = 0.012; % probability technology line shuts down

ggamma = .6; % match to labor share of GDP
epsilon = 2; % from Acemoglu, Akcigit, Bloom and Kerr (2013) - pg 21

rr = 0.05; % interest rate
bbeta = 1/(1+rr); % intertemporal discount factor
LL = 1; % production worker
PPi = ggamma/(1-bbeta); % constant part of the price of innovation

nu = 2.8; % chosen so that nrof firms makes avg growth rate is approximately 2%
nrofinv = exp(6); % number of inventors, chosen to match initial number of patents
nroffirms = round(nu*nrofinv); 
g1 = 0.066; % initial growth rate of patent numbers
g2 = 0.02; % final growth rate of patents

% define parameters for first 180 years
T_past = 180; % 1836 - 2016
params = v2struct(etaH, etaM, etaL, ttau, llambda, kkappa, xxi, zeta, ggamma, epsilon, rr, nu, g1, g2);

% run simulation from 1836 to 2016
result_before = simulate_path_zeta(T_past, params, nrofinv, 'before');

% sotore simulation results
quality_before = result_before.quality;
summat_before = result_before.summat;
state_before = result_before.state;
Mmat_before = result_before.Mmat;
Welfare_before = result_before.Welfare;
nrofinv_before = result_before.nrofinv;

% define parameters for counterfactuals
T_forward = 50;
subsidy_vals = linspace(.05,1,20);
subs_len = length(subsidy_vals);
subsidy = [0, 0; subsidy_vals', zeros(subs_len,1) ; zeros(subs_len,1), subsidy_vals']; 
    % subsidy for new tech, subsidy for new combo

% initialize matrices to store results
Welfare_after = zeros(1, size(subsidy,1));

% run simulation for counterfactual
for i = 1:size(subsidy,1)
    result_after = simulate_path_zeta(T_forward, params, nrofinv_before, 'after', subsidy(i,:), ...
    	quality_before, Mmat_before, summat_before, state_before);
    Welfare_after(i) = result_after.Welfare; 
end

Welfare = [repmat(Welfare_before, 1,  size(subsidy,1)); Welfare_after];

Welfare0 = Welfare_after(1)*ones(1,subs_len);
WelfareNT = Welfare_after(2:subs_len+1);
WelfareNC = Welfare_after(subs_len+2:end);

WChangeNT = -(WelfareNT - Welfare0) ./ Welfare0;
WChangeNC = -(WelfareNC - Welfare0) ./ Welfare0;

welfare_gains_pp = figure(1);
plot(subsidy_vals, WChangeNT, '--b', subsidy_vals, WChangeNC, 'r')
legend('New technologies', 'New combination', 'location', 'Northwest')
title('Welfare gains (%) associated with subsidy value')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_pp.png')

welfare_gains_NT = figure(2);
plot(subsidy_vals, WelfareNT, '--b', subsidy_vals, Welfare0, 'r')
legend('With subsidy', 'Baseline (so subsidy)', 'location', 'Northwest')
title('Welfare gains associated with subsidy for new technologies')
xlabel('Subsidy value')
ylabel('Welfare value')
saveas(gcf, 'tex_files/figures/welfare_NT.png')

welfare_gains_NC = figure(3);
plot(subsidy_vals, WelfareNC, '--b', subsidy_vals, Welfare0, 'r')
legend('With subsidy', 'Baseline (no subsidy)', 'location', 'Northwest')
title('Welfare gains associated with subsidy for new combinaitons')
xlabel('Subsidy value')
ylabel('Welfare value')
saveas(gcf, 'tex_files/figures/welfare_NC.png')






