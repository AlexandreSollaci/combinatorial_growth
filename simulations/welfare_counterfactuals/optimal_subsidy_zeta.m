clear
close all
clc

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/welfare_counterfactuals/')

etaH = 0.15; % NT step size
etaM = 0.06; % NC step size
etaL = 0; % Refinement step size
ttau = 1000; % shape parameter for ideas distribution
llambda = 1; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkappa = 2; % and shape parameter kappa
xxi = 200;   % 1/xi is the fraction of feasible combinations
zeta = 0.012; % probability technology line shuts down
ggamma = .6; % match to labor share of GDP
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

seed = 10;

% define parameters for iterations
Tmax = 230; % total number of periods
T_subs = 180; % start subsidy period in 2016
seed = 10;
params = v2struct(etaH, etaM, etaL, ttau, llambda, kkappa, xxi, zeta, ggamma, epsilon, rr, nu, nrofinv, g1, g2);

% Run simulation with no subsidy
subsidy_plain = [0,0];
result_plain = simulate_path_zeta(Tmax, T_subs, params, subsidy_plain, seed);

%% store simulation results
%Mmat_plain = result_plain.Mmat;
%summat_plain = result_plain.summat;
% GDP_plain = result_plain.GDP;
% Consumption_plain = result_plain.Consumption;
% inventor_cost_plain = result_plain.inventor_cost;
% firm_cost_plain = result_plain.firm_cost;
% policy_cost_plain = result_plain.policy_cost;
% nrofpatents_plain = result_plain.nrofpatents;
% quality_plain = result_plain.quality;
% Growth_plain = result_plain.Growth;
Welfare_plain = result_plain.Welfare;
Welfare_log_plain = result_plain.Welfare_log;

%% Run simulation with subsidy
subsidy_vals = [0.01, 0.02, 0.03, 0.04, linspace(.05,1,20)];
subs_len = length(subsidy_vals);
subsidy_subs = [subsidy_vals', zeros(subs_len,1)]; % ; zeros(subs_len,1), subsidy_vals'];

% initialize matrices
% GDP_subs = zeros(Tmax, subs_len);
% Consumption_subs = zeros(Tmax, subs_len);
% inventor_cost_subs = zeros(Tmax, subs_len);
% firm_cost_subs = zeros(Tmax, subs_len);
% policy_cost_subs = zeros(Tmax, subs_len);
% nrofpatents_subs = zeros(Tmax, subs_len);
% %quality_subs = result_subs.quality;
% Growth_subs = zeros(Tmax, subs_len);
Welfare_subs = zeros(subs_len,1);
Welfare_log_subs = zeros(subs_len,1);

for i = 1:subs_len
	result_subs = simulate_path_zeta(Tmax, T_subs, params, subsidy_subs(i,:), seed);

	% store simulation results
	% GDP_subs(:,i) = result_subs.GDP;
	% Consumption_subs(:,i) = result_subs.Consumption;
	% inventor_cost_subs(:,i) = result_subs.inventor_cost;
	% firm_cost_subs(:,i) = result_subs.firm_cost;
	% policy_cost_subs(:,i) = result_subs.policy_cost;
	% nrofpatents_subs(:,i) = result_subs.nrofpatents;
	% %quality_subs = result_subs.quality;
	% Growth_subs(:,i) = result_subs.Growth;
	Welfare_subs(i) = result_subs.Welfare;
	Welfare_log_subs(i) = result_subs.Welfare_log;
end

periods = 1836 + linspace(1,Tmax,Tmax);
periods_subs = 1836 + T_subs + linspace(1,Tmax-T_subs,Tmax-T_subs)';

WChangeNT = sign(Welfare_subs).*(Welfare_subs - Welfare_plain)/Welfare_plain;
WlogChangeNT = sign(Welfare_log_subs).*(Welfare_log_subs - Welfare_log_plain)/Welfare_log_plain;

welfare_gains_pp = figure(1);
plot(subsidy_vals, WChangeNT)
%legend('New technologies', 'New combination', 'location', 'Northeast')
title('Welfare gains (%) associated with subsidy value (CRRA)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_pp.png')

welfare_gains_pp = figure(1);
plot(subsidy_vals, WlogChangeNT)
%legend('New technologies', 'New combination', 'location', 'Northeast')
title('Welfare gains (%) associated with subsidy value (Log)')
xlabel('Subsidy value')
ylabel('Welfare gain')
saveas(gcf, 'tex_files/figures/welfare_log_pp.png')



























% % run simulation from 1836 to 2016
% result_before = simulate_path_zeta(T_past, params, nrofinv, 'before');

% % sotore simulation results
% quality_before = result_before.quality;
% summat_before = result_before.summat;
% state_before = result_before.state;
% Mmat_before = result_before.Mmat;
% Welfare_before = result_before.Welfare;
% nrofinv_before = result_before.nrofinv;

% % define parameters for counterfactuals
% T_forward = 50;
% subsidy_vals = [0.01, 0.02, 0.03, 0.04, linspace(.05,1,20)];
% subs_len = length(subsidy_vals);
% subsidy = [0, 0; subsidy_vals', zeros(subs_len,1) ; zeros(subs_len,1), subsidy_vals']; 
%     % subsidy for new tech, subsidy for new combo

% % initialize matrices to store results
% Welfare_after = zeros(1, size(subsidy,1));

% % run simulation for counterfactual
% for i = 1:size(subsidy,1)
%     result_after = simulate_path_zeta(T_forward, params, nrofinv_before, 'after', subsidy(i,:), ...
%     	quality_before, Mmat_before, summat_before, state_before);
%     Welfare_after(i) = result_after.Welfare; 
% end

% Welfare = [repmat(Welfare_before, 1,  size(subsidy,1)); Welfare_after];

% Welfare0 = Welfare_after(1)*ones(1,subs_len);
% WelfareNT = Welfare_after(2:subs_len+1);
% WelfareNC = Welfare_after(subs_len+2:end);

% WChangeNT = -(WelfareNT - Welfare0) ./ Welfare0;
% WChangeNC = -(WelfareNC - Welfare0) ./ Welfare0;

% welfare_gains_pp = figure(1);
% plot(subsidy_vals, WChangeNT, '--b', subsidy_vals, WChangeNC, 'r')
% legend('New technologies', 'New combination', 'location', 'Northeast')
% title('Welfare gains (%) associated with subsidy value')
% xlabel('Subsidy value')
% ylabel('Welfare gain')
% saveas(gcf, 'tex_files/figures/welfare_pp.png')

% welfare_gains_NT = figure(2);
% plot(subsidy_vals, WelfareNT, '--b', subsidy_vals, Welfare0, 'r')
% legend('With subsidy', 'Baseline (no subsidy)', 'location', 'Northeast')
% title('Welfare gains associated with subsidy for new technologies')
% xlabel('Subsidy value')
% ylabel('Welfare value')
% saveas(gcf, 'tex_files/figures/welfare_NT.png')

% welfare_gains_NC = figure(3);
% plot(subsidy_vals, WelfareNC, '--b', subsidy_vals, Welfare0, 'r')
% legend('With subsidy', 'Baseline (no subsidy)', 'location', 'Southwest')
% title('Welfare gains associated with subsidy for new combinaitons')
% xlabel('Subsidy value')
% ylabel('Welfare value')
% saveas(gcf, 'tex_files/figures/welfare_NC.png')






