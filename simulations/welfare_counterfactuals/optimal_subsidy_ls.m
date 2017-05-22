clear
close all
clc

feature accel on

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/welfare_counterfactuals/')

etaH = 0.15; % NT step size
etaM = 0.015; % NC step size
etaL = 0; % Refinement step size
ttau = 500; % shape parameter for ideas distribution
phi = 1.075;
llambda = 2; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkappa = 2;  % and shape parameter kappa
xxi = 5;   % 1/xi is the fraction of feasible combinations
zeta = 0.01; % probability technology line shuts down
ggamma = 0.6; % match to labor share of GDP
epsilon = 2; % from Acemoglu, Akcigit, Bloom and Kerr (2013) - pg 21
rr = 0.05; % interest rate
bbeta = 1/(1+rr); % intertemporal discount factor
LL = 1; % production worker
PPi = ggamma/(1-bbeta); % constant part of the price of innovation
nu = 2; % chosen so that nrof firms makes avg growth rate is approximately 2%
nrofinv = exp(6); % initial number of inventors, chosen to match initial number of patents
nroffirms = round(nu*nrofinv); 
g1 = 0.066; % initial growth rate of patent numbers
g2 = 0.02; % final growth rate of patents

% define parameters for iterations
Tmax = 50; % total number of periods
T_subs = 30; % start subsidy period in 2016
seed = 1702;
params = v2struct(etaH, etaM, etaL, ttau, phi, llambda, kkappa, xxi, zeta, ggamma, epsilon, rr, nu, nrofinv, g1, g2);

%% Run simulation with no subsidy
subsidy_plain = [0,0];
result_plain = simulate_path_ls(Tmax, T_subs, params, subsidy_plain, seed);

Welfare_plain = result_plain.Welfare;
Welfare_log_plain = result_plain.Welfare_log;

%% Run simulation with subsidy
subsidy_vals = linspace(0,etaH,5);
subs_len = length(subsidy_vals);

% new technologies
subsidy_substech = [subsidy_vals', zeros(subs_len,1)]; 
Welfare_substech = zeros(subs_len,1);
Welfare_log_substech = zeros(subs_len,1);

for i = 1:subs_len
	result_subs = simulate_path_ls(Tmax, T_subs, params, subsidy_substech(i,:), seed);

	Welfare_substech(i) = result_subs.Welfare;
	Welfare_log_substech(i) = result_subs.Welfare_log;
end

% new combinations
subsidy_subscomb = [zeros(subs_len,1), subsidy_vals']; 
Welfare_subscomb = zeros(subs_len,1);
Welfare_log_subscomb = zeros(subs_len,1);

for i = 1:subs_len
	result_subs = simulate_path_ls(Tmax, T_subs, params, subsidy_subscomb(i,:), seed);

	Welfare_subscomb(i) = result_subs.Welfare;
	Welfare_log_subscomb(i) = result_subs.Welfare_log;
end

% both
subsidy_subsboth = [subsidy_vals', subsidy_vals']; 
Welfare_subsboth = zeros(subs_len,1);
Welfare_log_subsboth = zeros(subs_len,1);

for i = 1:subs_len
	result_subs = simulate_path_ls(Tmax, T_subs, params, subsidy_subsboth(i,:), seed);

	Welfare_subsboth(i) = result_subs.Welfare;
	Welfare_log_subsboth(i) = result_subs.Welfare_log;
end

periods = 1836 + linspace(1,Tmax,Tmax);
periods_subs = 1836 + T_subs + linspace(1,Tmax-T_subs,Tmax-T_subs)';

WChangeNT = sign(Welfare_substech).*(Welfare_substech - Welfare_plain)/Welfare_plain;
WlogChangeNT = sign(Welfare_log_substech).*(Welfare_log_substech - Welfare_log_plain)/Welfare_log_plain;

WChangeNC = sign(Welfare_subscomb).*(Welfare_subscomb - Welfare_plain)/Welfare_plain;
WlogChangeNC = sign(Welfare_log_subscomb).*(Welfare_log_subscomb - Welfare_log_plain)/Welfare_log_plain;

WChangeTC = sign(Welfare_subsboth).*(Welfare_subsboth - Welfare_plain)/Welfare_plain;
WlogChangeTC = sign(Welfare_log_subsboth).*(Welfare_log_subsboth - Welfare_log_plain)/Welfare_log_plain;

welfare_gains_pp = figure(1);
plot(subsidy_vals, WChangeNT, '--r', subsidy_vals, WChangeNC, '-.b', subsidy_vals, WChangeTC, 'k')
legend('New technologies', 'New combination', 'Both', 'location', 'Northwest')
title('Welfare gains (%) associated with subsidy value (CRRA)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_ls_3ways.png')

welfare_log_gains_pp = figure(2);
plot(subsidy_vals, WlogChangeNT, '--r', subsidy_vals, WlogChangeNC, '-.b', subsidy_vals, WlogChangeTC, 'k')
legend('New technologies', 'New combination', 'Both', 'location', 'Northwest')
title('Welfare gains (%) associated with subsidy value (Log)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_log_ls_3ways.png')

welfare_gains_pp_NT = figure(3);
plot(subsidy_vals, WChangeNT, 'k')
%legend('New technologies', 'New combination', 'Both', 'location', 'Northwest')
title('Welfare gains (%) associated with new tech subsidy (CRRA)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_ls_NT.png')

welfare_log_gains_pp_NT = figure(4);
plot(subsidy_vals, WlogChangeNT, 'k')
%legend('New technologies', 'New combination', 'Both', 'location', 'Northwest')
title('Welfare gains (%) associated with new tech subsidy (Log)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_log_ls_NT.png')
