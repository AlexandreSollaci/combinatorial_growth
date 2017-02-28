clear
close all
clc

%cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/welfare_counterfactuals/')

% parameters
etaH = 0.15; % NT step size
etaM = 0.075; % NC step size
etaL = 0; % Refinement step size
ttau = 500; % shape parameter for ideas distribution
phi = 1.05;
llambda = 1.5; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkappa = 2;  % and shape parameter kappa
xxi = 175;   % 1/xi is the fraction of feasible combinations
zeta = 0.01; % probability technology line shuts down
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

% define parameters for first 180 years
Tmax = 280; % 1836 - 2016
T_subs = 180; % start subsidy period in 2016
seed = 10;
params = v2struct(etaH, etaM, etaL, ttau, phi, llambda, kkappa, xxi, zeta, ggamma, epsilon, rr, nu, nrofinv, g1, g2);

% Run simulation with no subsidy
subsidy_plain = [0,0];
result_plain = simulate_path_zeta(Tmax, T_subs, params, subsidy_plain, seed);

% store simulation results
Mmat_plain = result_plain.Mmat;
summat_plain = result_plain.summat;
GDP_plain = result_plain.GDP;
Consumption_plain = result_plain.Consumption;
inventor_cost_plain = result_plain.inventor_cost;
firm_cost_plain = result_plain.firm_cost;
policy_cost_plain = result_plain.policy_cost;
nrofpatents_plain = result_plain.nrofpatents;
quality_plain = result_plain.quality;
Growth_plain = result_plain.Growth;
Welfare_plain = result_plain.Welfare;
Welfare_log_plain = result_plain.Welfare_log;

% Run simulation with subsidy
subsidy_subs = [0.2,0];
result_subs = simulate_path_zeta(Tmax, T_subs, params, subsidy_subs, seed);

% store simulation results
Mmat_subs = result_subs.Mmat;
summat_subs = result_subs.summat;
GDP_subs = result_subs.GDP;
Consumption_subs = result_subs.Consumption;
inventor_cost_subs = result_subs.inventor_cost;
firm_cost_subs = result_subs.firm_cost;
policy_cost_subs = result_subs.policy_cost;
nrofpatents_subs = result_subs.nrofpatents;
quality_subs = result_subs.quality;
Growth_subs = result_subs.Growth;
Welfare_subs = result_subs.Welfare;
Welfare_log_subs = result_subs.Welfare_log;

% number of new technologies produced each year
numtech_plain = nrofpatents_plain.*summat_plain(2:end,2);
numtech_subs = nrofpatents_subs.*summat_subs(2:end,2);

periods = 1836 + linspace(1,Tmax,Tmax);
periods_subs = 1836 + T_subs + linspace(1,Tmax-T_subs,Tmax-T_subs)';

aggregates = figure(1);
subplot(3,2,1)
plot(periods, log(Consumption_plain), 'r', periods, log(Consumption_subs), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Consumption')
xlabel('Year'); ylabel('Consumption')

subplot(3,2,2)
plot(periods, log(GDP_plain), 'r', periods, log(GDP_subs), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of GDP')
xlabel('Year'); ylabel('GDP')

subplot(3,2,3)
plot(periods, log(firm_cost_plain), 'r', periods, log(firm_cost_subs), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Firm costs')
xlabel('Year'); ylabel('Firms` Costs')

subplot(3,2,4)
plot(periods, log(inventor_cost_plain), 'r', periods, log(inventor_cost_subs), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Inventor costs')
xlabel('Year') ;ylabel('Inventors` Costs')

subplot(3,2,5)
plot(periods, log(policy_cost_plain), 'r', periods, log(policy_cost_subs), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Policy costs')
xlabel('Year') ;ylabel('Subsidy Costs')

subplot(3,2,6)
plot(periods, Growth_plain, 'r', periods, Growth_subs, '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northeast')
title('Growth rate of the economy')
xlabel('Year'); ylabel('Growth rate')
xlim([1836, 1836+Tmax])
saveas(gcf, 'tex_files/figures/aggregates.png')

aggregates2016 = figure(2);
subplot(3,2,1)
plot(periods_subs, log(Consumption_plain(T_subs+1:end)), 'r', periods_subs, log(Consumption_subs(T_subs+1:end)), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Consumption')
xlabel('Year'); ylabel('Consumption')
xlim([1836+T_subs, 1836+Tmax])

subplot(3,2,2)
plot(periods_subs, log(GDP_plain(T_subs+1:end)), 'r', periods_subs, log(GDP_subs(T_subs+1:end)), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of GDP')
xlabel('Year'); ylabel('GDP')
xlim([1836+T_subs, 1836+Tmax])

subplot(3,2,3)
plot(periods_subs, log(firm_cost_plain(T_subs+1:end)), 'r', periods_subs, log(firm_cost_subs(T_subs+1:end)), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Firm costs')
xlabel('Year'); ylabel('Firms` Costs')
xlim([1836+T_subs, 1836+Tmax])

subplot(3,2,4)
plot(periods_subs, log(inventor_cost_plain(T_subs+1:end)), 'r', periods_subs, log(inventor_cost_subs(T_subs+1:end)), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Inventor costs')
xlabel('Year') ;ylabel('Inventors` Costs')
xlim([1836+T_subs, 1836+Tmax])

subplot(3,2,5)
plot(periods_subs, log(policy_cost_plain(T_subs+1:end)), 'r', periods_subs, log(policy_cost_subs(T_subs+1:end)), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Policy costs')
xlabel('Year') ;ylabel('Subsidy Costs')
xlim([1836+T_subs, 1836+Tmax])

subplot(3,2,6)
plot(periods_subs, Growth_plain(T_subs+1:end), 'r', periods_subs, Growth_subs(T_subs+1:end), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Growth rate of the economy')
xlabel('Year'); ylabel('Growth rate')
xlim([1836+T_subs, 1836+Tmax])
saveas(gcf, 'tex_files/figures/aggregates2016.png')

aggregates2016_diff = figure(3);
subplot(2,2,1)
plot(periods_subs, log(Consumption_subs(T_subs+1:end)) - log(Consumption_plain(T_subs+1:end)))
title('Change in log of Consumption')
xlabel('Year'); ylabel('Consumption')
xlim([1836+T_subs, 1836+Tmax])

subplot(2,2,2)
plot(periods_subs, log(GDP_subs(T_subs+1:end)) - log(GDP_plain(T_subs+1:end)))
title('Change in log of GDP')
xlabel('Year'); ylabel('GDP')
xlim([1836+T_subs, 1836+Tmax])

subplot(2,2,3)
plot(periods_subs, log(firm_cost_subs(T_subs+1:end)) - log(firm_cost_plain(T_subs+1:end)))
title('Change in log of Firm costs')
xlabel('Year'); ylabel('Firms` Costs')
xlim([1836+T_subs, 1836+Tmax])

subplot(2,2,4)
plot(periods_subs, log(inventor_cost_subs(T_subs+1:end)) - log(inventor_cost_plain(T_subs+1:end)))
title('Change in log of Inventor costs')
xlabel('Year') ;ylabel('Inventors` Costs')
xlim([1836+T_subs, 1836+Tmax])
saveas(gcf, 'tex_files/figures/aggregates2016_diff.png')

patents = figure(4);
subplot(2,2,1)
plot(summat_plain(:,1), summat_plain(:,2), 'r', summat_subs(:,1), summat_subs(:,2), '--b')
legend('No subsidy', 'Subsdidy')
title('New technology evolution')
xlabel('Year'); ylabel('Fraction of new techs')
xlim([1836, 1836+Tmax])

subplot(2,2,2)
plot(summat_plain(:,1), summat_plain(:,3), 'r', summat_subs(:,1), summat_subs(:,3), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Southeast')
title('New combination evolution')
xlabel('Year'); ylabel('Fraction of new combos')
xlim([1836, 1836+Tmax])

subplot(2,2,3)
plot(summat_plain(:,1), summat_plain(:,4), 'r', summat_subs(:,1), summat_subs(:,4), '--b')
legend('No subsidy', 'Subsdidy')
title('Reuse/refinement evolution')
xlabel('Year'); ylabel('Fraction of reuse')
xlim([1836, 1836+Tmax])

subplot(2,2,4)
plot(periods, numtech_plain, 'r', periods, numtech_subs, '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Number of technologies produced each year')
xlabel('Year'); ylabel('Number of patents')
xlim([1836, 1836+Tmax])
saveas(gcf, 'tex_files/figures/patents.png')


patents = figure(5);
subplot(2,2,1)
plot(summat_plain(T_subs+1:end,1), summat_plain(T_subs+1:end,2), 'r', summat_subs(T_subs+1:end,1), summat_subs(T_subs+1:end,2), '--b')
legend('No subsidy', 'Subsdidy')
title('New technology evolution')
xlabel('Year'); ylabel('Fraction of new techs')
xlim([1836+T_subs, 1836+Tmax])
ylim([0,.02])

subplot(2,2,2)
plot(summat_plain(T_subs+1:end,1), summat_plain(T_subs+1:end,3), 'r', summat_subs(T_subs+1:end,1), summat_subs(T_subs+1:end,3), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Southwest')
title('New combination evolution')
xlabel('Year'); ylabel('Fraction of new combos')
xlim([1836+T_subs, 1836+Tmax])

subplot(2,2,3)
plot(summat_plain(T_subs+1:end,1), summat_plain(T_subs+1:end,4), 'r', summat_subs(T_subs+1:end,1), summat_subs(T_subs+1:end,4), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Reuse/refinement evolution')
xlabel('Year'); ylabel('Fraction of reuse')
xlim([1836+T_subs, 1836+Tmax])

subplot(2,2,4)
plot(periods_subs, numtech_plain(T_subs+1:end), 'r', periods_subs, numtech_subs(T_subs+1:end), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Number of technologies produced each year')
xlabel('Year'); ylabel('Number of patents')
xlim([1836+T_subs, 1836+Tmax])
saveas(gcf, 'tex_files/figures/patents2016.png')


Welfare_gain = sign(Welfare_plain)*(Welfare_subs - Welfare_plain) / Welfare_plain;
Welfare_log_gain = sign(Welfare_log_plain)*(Welfare_log_subs - Welfare_log_plain) / Welfare_log_plain;

display(['Welfare gain with CRRA utility is ' num2str(Welfare_gain)])
display(['Welfare gain with log utlity is ' num2str(Welfare_log_gain)])















% % define parameters for counterfactuals
% T_forward = 50;
% subsidy = [0, 0; -10, 0]; % sibsidy for new tech, subsidy for new combo

% % initialize matrices to store results
% Consumption_after = zeros(T_forward, size(subsidy,1));
% GDP_after = zeros(T_forward, size(subsidy,1));
% inventor_cost_after = zeros(T_forward, size(subsidy,1));
% firm_cost_after = zeros(T_forward, size(subsidy,1));
% policy_cost_after = zeros(T_forward, size(subsidy,1));
% Discounted_utility_after = zeros(T_forward, size(subsidy,1));
% Growth_after = zeros(T_forward, size(subsidy,1));
% Welfare_after = zeros(1, size(subsidy,1));
% summat_after = zeros(T_forward + 1, 4,  size(subsidy,1));

% numtech = zeros(T_past + T_forward, size(subsidy,1));
% % run simulation for counterfactual
% seed = 10;
% for i = 1:size(subsidy,1)
%     result_after = simulate_path_zeta(T_forward, params, nrofinv_before, 'after', subsidy(i,:), ...
%         quality_before, Mmat_before, summat_before, state_before, seed);

%     Consumption_after(:,i) = result_after.Consumption;
%     GDP_after(:,i) = result_after.GDP;
%     inventor_cost_after(:,i) = result_after.inventor_cost;
%     firm_cost_after(:,i) = result_after.firm_cost; 
%     policy_cost_after(:,i) = result_after.policy_cost;
%     Discounted_utility_after(:,i) = result_after.Discounted_utility;
%     Growth_after(:,i) = result_after.Growth; 
%     Welfare_after(i) = result_after.Welfare; 
%     summat_after(:,:,i) = result_after.summat; % [year, new tech, new combo, reuse]

%     Mmat_after = result_after.Mmat;

%     for t = 1:T_past+T_forward
%         index = find(Mmat_after(:,4) == 1836 + t, 1, 'last');
%         if t <= T_past
%             numtech(t, :) = repmat(length(Mmat_after(1:index, 1)), 1, size(subsidy,1));
%         else
%             numtech(t,i) = length(Mmat_after(1:index, 1));
%         end
%     end

% end

% Consumption = [repmat(Consumption_before, 1,  size(subsidy,1)); Consumption_after];
% GDP = [repmat(GDP_before, 1,  size(subsidy,1)); GDP_after];
% inventor_cost = [repmat(inventor_cost_before, 1,  size(subsidy,1)); inventor_cost_after];
% firm_cost = [repmat(firm_cost_before, 1,  size(subsidy,1)); firm_cost_after];
% policy_cost = [repmat(policy_cost_before, 1,  size(subsidy,1)); policy_cost_after];
% Discounted_utility = [repmat(Discounted_utility_before, 1,  size(subsidy,1)); Discounted_utility_after];
% Growth = [repmat(Growth_before, 1,  size(subsidy,1)); Growth_after];
% Welfare = [repmat(Welfare_before, 1,  size(subsidy,1)); Welfare_after];
% summat = [repmat(summat_before(:,:), 1, 1, size(subsidy, 1)) ; summat_after(2:end,:,:)];

% Tmax = T_past + T_forward;
% periods = 1836 + linspace(1,Tmax,Tmax);
% periods_after = 1836 + T_past + linspace(1, T_forward, T_forward);

% aggregates = figure(1);
% subplot(3,2,1)
% plot(periods, log(Consumption(:,1)), 'r', periods, log(Consumption(:,2)), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Consumption')
% xlabel('Year'); ylabel('Consumption')

% subplot(3,2,2)
% plot(periods, log(GDP(:,1)), 'r', periods, log(GDP(:,2)), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of GDP')
% xlabel('Year'); ylabel('GDP')

% subplot(3,2,3)
% plot(periods, log(firm_cost(:,1)), 'r', periods, log(firm_cost(:,2)), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Firm costs')
% xlabel('Year'); ylabel('Firms` Costs')

% subplot(3,2,4)
% plot(periods, log(inventor_cost(:,1)), 'r', periods, log(inventor_cost(:,2)), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Inventor costs')
% xlabel('Year') ;ylabel('Inventors` Costs')

% subplot(3,2,5)
% plot(periods, log(policy_cost(:,1)), 'r', periods, log(policy_cost(:,2)), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Policy costs')
% xlabel('Year') ;ylabel('Subsidy Costs')

% subplot(3,2,6)
% plot(periods, Growth(:,1), 'r', periods, Growth(:,2), '--b')
% legend('No subsidy', 'Subsdidy', 'location', 'Northeast')
% title('Growth rate of the economy')
% xlabel('Year'); ylabel('Growth rate')
% xlim([1836, 1836+Tmax])
% saveas(gcf, 'tex_files/figures/aggregates.png')

% aggregates2016 = figure(2);
% subplot(3,2,1)
% plot(periods_after, log(Consumption_after(:,1)), 'r', periods_after, log(Consumption_after(:,2)), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Consumption, starting in 2016')
% xlabel('Year'); ylabel('Consumption')

% subplot(3,2,2)
% plot(periods_after, log(GDP_after(:,1)), 'r', periods_after, log(GDP_after(:,2)), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of GDP, starting in 2016')
% xlabel('Year'); ylabel('GDP')

% subplot(3,2,3)
% plot(periods_after, log(firm_cost_after(:,1)), 'r', periods_after, log(firm_cost_after(:,2)), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Firm costs, starting in 2016')
% xlabel('Year'); ylabel('Firms` Costs')

% subplot(3,2,4)
% plot(periods_after, log(inventor_cost_after(:,1)), 'r', periods_after, log(inventor_cost_after(:,2)), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Inventor costs, starting in 2016')
% xlabel('Year') ;ylabel('Inventors` Costs')

% subplot(3,2,5)
% plot(periods_after, log(policy_cost_after(:,1)), 'r', periods_after, log(policy_cost_after(:,2)), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('log of Policy costs, starting in 2016')
% xlabel('Year') ;ylabel('Subsidy Costs')

% subplot(3,2,6)
% plot(periods_after, Growth_after(:,1), 'r', periods_after, Growth_after(:,2), '--b')
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('Growth rate of the economy, starting in 2016')
% xlim([1836+T_past, 1836+Tmax])
% saveas(gcf, 'tex_files/figures/aggregates2016.png')

% aggregates_diff = figure(3);
% subplot(3,1,1)
% plot(periods_after, log(Consumption_after(:,2)) - log(Consumption_after(:,1)))
% xlim([1836+T_past, 1836+Tmax])
% title('log diff of Consumption, with and without subsidy')
% xlabel('Year'); ylabel('log difference of Consumption')

% subplot(3,1,2)
% plot(periods_after, log(GDP_after(:,2)) - log(GDP_after(:,1)))
% xlim([1836+T_past, 1836+Tmax])
% title('log diff of GDP, with and without subsidy')
% xlabel('Year'); ylabel('log difference of GDP')

% subplot(3,1,3)
% plot(periods_after, log(firm_cost_after(:,2)) - log(firm_cost_after(:,1)))
% xlim([1836+T_past, 1836+Tmax])
% title('log diff of Firms costs, with and without subsidy')
% xlabel('Year'); ylabel('log difference of Firms costs')
% saveas(gcf, 'tex_files/figures/aggregates_diff.png')

% newtechcreation = figure(4);
% subplot(2,1,1)
% plot(periods(2:end), diff(numtech(:,1)), 'r', periods(2:end), diff(numtech(:,2)), '--b')
% legend('No subsidy', 'Subsdidy')
% title('Number of technologies created')
% xlabel('Year'); ylabel('Number of new techs')
% xlim([1836, 1836+Tmax])

% subplot(2,1,2)
% plot(periods_after(1:end), diff(numtech(T_past:end,1)), 'r', periods_after(1:end), diff(numtech(T_past:end,2)), '--b')
% legend('No subsidy', 'Subsdidy')
% title('Number of technologies created, starting in 2016')
% xlabel('Year'); ylabel('Number of new techs')
% xlim([1836+T_past, 1836+Tmax])
% saveas(gcf, 'tex_files/figures/tech_creation.png')

% patents = figure(5);
% subplot(2,2,1)
% plot(summat(:,1,1), summat(:,2,1), 'r', summat(:,1,2), summat(:,2,2), '--b')
% legend('No subsidy', 'Subsdidy')
% title('New technology evolution')
% xlabel('Year'); ylabel('Fraction of new techs')
% xlim([1836, 1836+Tmax])

% subplot(2,2,2)
% plot(summat(:,1,1), summat(:,3,1), 'r', summat(:,1,2), summat(:,3,2), '--b')
% legend('No subsidy', 'Subsdidy', 'location', 'Southeast')
% title('New combination evolution')
% xlabel('Year'); ylabel('Fraction of new combos')
% xlim([1836, 1836+Tmax])

% subplot(2,2,3)
% plot(summat(:,1,1), summat(:,4,1), 'r', summat(:,1,2), summat(:,4,2), '--b')
% legend('No subsidy', 'Subsdidy')
% title('Reuse/refinement evolution')
% xlabel('Year'); ylabel('Fraction of reuse')
% xlim([1836, 1836+Tmax])

% subplot(2,2,4)
% plot(periods, numtech(:,1), 'r', periods, numtech(:,2), '--b')
% xlim([1836, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('Number of existing technologies')
% xlabel('Year') ;ylabel('Number of product lines')
% saveas(gcf, 'tex_files/figures/patents.png')

% patents2016 = figure(6);
% subplot(2,2,1)
% plot(summat_after(:,1,1), summat_after(:,2,1), 'r', summat_after(:,1,2), summat_after(:,2,2), '--b')
% legend('No subsidy', 'Subsdidy')
% title('New technology evolution, starting in 2016')
% xlabel('Year'); ylabel('Fraction of new techs')
% xlim([1836+T_past, 1836+Tmax])

% subplot(2,2,2)
% plot(summat_after(:,1,1), summat_after(:,3,1), 'r', summat_after(:,1,2), summat_after(:,3,2), '--b')
% legend('No subsidy', 'Subsdidy', 'location', 'Southeast')
% title('New combination evolution, starting in 2016')
% xlabel('Year'); ylabel('Fraction of new combos')
% xlim([1836+T_past, 1836+Tmax])

% subplot(2,2,3)
% plot(summat_after(:,1,1), summat_after(:,4,1), 'r', summat_after(:,1,2), summat_after(:,4,2), '--b')
% legend('No subsidy', 'Subsdidy')
% title('Reuse/refinement evolution, starting in 2016')
% xlabel('Year'); ylabel('Fraction of reuse')
% xlim([1836+T_past, 1836+Tmax])

% subplot(2,2,4)
% plot(periods_after, numtech(T_past+1:end,1), 'r', periods_after, numtech(T_past+1:end,2), '--b')
% xlim([1836+T_past, 1836+Tmax])
% legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
% title('Number of existing technologies, starting in 2016')
% xlabel('Year') ;ylabel('Number of product lines')
% saveas(gcf, 'tex_files/figures/patents2016.png')

% T = table;
% T.Welfare_no_sub = Welfare_after(1);
% T.Welfare_sub = Welfare_after(2);
% T.Welfare_increase = sign(Welfare_after(1))*( Welfare_after(2) - Welfare_after(1) ) / Welfare_after(1);
% T.Numtech_no_sub = numtech(end,1);
% T.Numtech_sub = numtech(end,2);
% T.Numtech_increase = (numtech(end,2) - numtech(end,1)) / numtech(end,1);

% T2 = table;
% T2.Average_growth_no_sub = mean(Growth_after(:,1));
% T2.Average_growth_sub = mean(Growth_after(:,2));
% T2.Average_growth_increase = ( mean(Growth_after(:,2)) - mean(Growth_after(:,1)) )/ mean(Growth_after(:,1));

% display(T)
% display(T2)


% Summary_stats = [Welfare_after , sign(Welfare_after(1))*( Welfare_after(2) - Welfare_after(1) ) / Welfare_after(1) ; ...
%     numtech(end,:) , (numtech(end,2) - numtech(end,1)) / numtech(end,1); ...
%     mean(Growth_after) , ( mean(Growth_after(:,2)) - mean(Growth_after(:,1)) )/ mean(Growth_after(:,1)) ];

% rowLabels = {'Welfare'; '\# technologies'; 'Average GDP growth'};
% columnLabels = {'No Subsidy', 'Subsidy', '\% Change'};
% matrix2latex(Summary_stats, 'tex_files/counterfactual_table.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels)

