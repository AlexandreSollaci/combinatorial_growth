clear
close all
clc

%cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/consolidated/')

rng(10)

%% Parameters

etaH = 0.15;            % NT step size
etaM = 0.1;             % NC step size
etaL = 0;               % Refinement step size
tau = 500;              % shape parameter for ideas distribution
phi = 1.075;            % controls the mode of distribution of new ideas
lambda = 2;             % alpha (cost) has weibull distribution with scale parameter lambda ...
kappa = 2;              % and shape parameter kappa
xi = 5;                 % 1/xi is the fraction of feasible combinations
zeta = 0.01;            % probability technology line shuts down
gamma = 0.6;            % match to labor share of GDP
epsilon = 2;            % from Acemoglu, Akcigit, Bloom and Kerr (2013) - pg 21
beta = 1/(1.05);        % intertemporal discount factor
Pi = gamma/(1-beta);    % constant part of the price of innovation
nu = 3.1;               % ratio of firms to inventors; chosen so that nrof firms makes avg growth rate is approximately 2%
g1 = 0.066;             % initial growth rate of patent numbers
g2 = 0.02;              % final growth rate of patents

% store parameters
params = v2struct(etaH, etaM, etaL, tau, phi, lambda, kappa, xi, zeta, gamma, epsilon, beta, Pi, nu, g1, g2);

beg_inv = exp(6);       % number of inventors, chosen to match initial number of patents
beg_year = 1836;        % beginning year
Tbase = 10;            % number of periods to run baseline simulation
Tsubs = 10;             % number of periods to run subsidy simulation
Tmax = Tbase + Tsubs;

%% RUN THE MODEL

% Matrix to keep track of state of each product line
M0 = 10;                %initial number of technologies
% [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
Mmat0 = zeros(M0,5);
Mmat0(1:M0,1) = 1:M0;                % m ID
Mmat0(1:M0,2) = randi([1 3],M0,1);   % average tech pool in 1836 was 2.
Mmat0(1:M0,3) = zeros(M0,1);         % reuse/refinement so far
Mmat0(:,4) =  beg_year;              % birth year
Mmat0(1:M0,5) = ones(M0,1);          % technologies are born hot

% Matrix to keep track of fractions of each patent type
% [year, new tech, new comb, reuse]
patmat0 = zeros(Tbase,4); 
patmat1 = zeros(Tsubs,4); 

% state of product line: 1-new tech 2-recomb 3-reuse
state0 = zeros(1,M0);
state0(1:M0) = 1;

% Simulate Model

num_simul = 10;         % number of times model is simulated
patmat = zeros(Tmax+1, 4, num_simul);
growth = zeros(Tmax, num_simul);
consumption = zeros(Tmax, num_simul);
GDP = zeros(Tmax, num_simul);
firm_cost = zeros(Tmax, num_simul);
inventor_cost = zeros(Tmax, num_simul);
policy_cost = zeros(Tmax, num_simul);
nrofpatents = zeros(Tmax, num_simul);
welfare = zeros(2, num_simul);

subsidy_vals = [.05 .1 .2 .3];
subsidy = [subsidy_vals' zeros(length(subsidy_vals),1); ...
            zeros(length(subsidy_vals),1) subsidy_vals';...
            subsidy_vals' subsidy_vals']; % subsidy for [new tech, new comb]

patmat_subs = zeros(Tmax+1, 4, num_simul, length(subsidy));
growth_subs = zeros(Tmax, num_simul, length(subsidy));
consumption_subs = zeros(Tmax, num_simul, length(subsidy));
GDP_subs = zeros(Tmax, num_simul, length(subsidy));
firm_cost_subs = zeros(Tmax, num_simul, length(subsidy));
inventor_cost_subs = zeros(Tmax, num_simul, length(subsidy));
policy_cost_subs = zeros(Tmax, num_simul, length(subsidy));
nrofpatents_subs = zeros(Tmax, num_simul, length(subsidy));
welfare_subs = zeros(2, num_simul, length(subsidy));

for s = 1:num_simul
    seed1 = 1000*rand;
    model_plain = simulate_model(params, Tbase, beg_year, beg_inv, Mmat0, patmat0, state0, [0 0], seed1);

    Mmat = model_plain.Mmat;
    state = model_plain.state;
    quality = model_plain.quality;
    nrofinv = model_plain.nrofinv;

    seed2 = 1000*rand;
    
    model_subs_0 = simulate_model(params, Tsubs, beg_year+Tbase, nrofinv, Mmat, patmat1, state, [0 0], seed2);

    % In period 0, all ideas are new.
    patmat(:, :, s) = [[beg_year, 1, 0, 0]; model_plain.patmat; model_subs_0.patmat];
    nrofpatents(:,s) = [model_plain.nrofpatents; model_subs_0.nrofpatents];
    consumption(:,s) = [model_plain.consumption; model_subs_0.consumption];
    GDP(:,s) = [model_plain.GDP; model_subs_0.GDP];
    firm_cost(:,s) = [model_plain.firm_cost; model_subs_0.firm_cost];
    inventor_cost(:,s) = [model_plain.inventor_cost; model_subs_0.inventor_cost];
    policy_cost(:,s) = [model_plain.policy_cost; model_subs_0.policy_cost];
    growth(:,s) = [model_plain.growth; model_subs_0.growth];
    welfare(:,s) = [model_plain.welfare; model_subs_0.welfare];

    for v = 1:length(subsidy)

        model_subs_1 = simulate_model(params, Tsubs, beg_year+Tbase, nrofinv, Mmat, patmat1, state, subsidy(v,:), seed2);

        patmat_subs(:,:,s,v) = [[beg_year, 1, 0, 0]; model_plain.patmat; model_subs_1.patmat];
        nrofpatents_subs(:,s,v) = [model_plain.nrofpatents; model_subs_1.nrofpatents];
        consumption_subs(:,s,v) = [model_plain.consumption; model_subs_1.consumption];
        GDP_subs(:,s,v) = [model_plain.GDP; model_subs_1.GDP];
        firm_cost_subs(:,s,v) = [model_plain.firm_cost; model_subs_1.firm_cost];
        inventor_cost_subs(:,s,v) = [model_plain.inventor_cost; model_subs_1.inventor_cost];
        policy_cost_subs(:,s,v) = [model_plain.policy_cost; model_subs_1.policy_cost];
        growth_subs(:,s,v) = [model_plain.growth; model_subs_1.growth];
        welfare_subs(:,s,v) = [model_plain.welfare; model_subs_1.welfare];

    end
end

patent_shares = mean(patmat,3);
nrofpatents = mean(nrofpatents,2);
GDP = mean(GDP,2);
firm_cost = mean(firm_cost,2);
consumption = mean(consumption,2);
inventor_cost = mean(inventor_cost,2);
policy_cost = mean(policy_cost,2);
growth = mean(growth,2);
welfare = mean(welfare,2);

patent_shares_subs = mean(patmat_subs,3);
nrofpatents_subs = mean(nrofpatents_subs,2);
GDP_subs = mean(GDP_subs,2);
firm_cost_subs = mean(firm_cost_subs,2);
consumption_subs = mean(consumption_subs,2);
inventor_cost_subs = mean(inventor_cost_subs,2);
policy_cost_subs = mean(policy_cost_subs,2);
growth_subs = mean(growth_subs,2);
welfare_subs = mean(welfare_subs,2);

%% PLOTS

% Summary of Economy
periods = beg_year + linspace(1,Tmax,Tmax);
figure(1)
plot(periods, growth, '-k')
xlim([beg_year, beg_year+Tmax])
title('Rate of Growth')
saveas(gcf, 'figures/growth.png')

figure(2)
plot(periods, log(GDP), 'k', periods, log(consumption), '--k', periods, log(firm_cost), '-^k', periods, log(inventor_cost), '-ok')
legend('GDP', 'Consumption', 'Firms` Cost', 'Inventors` Cost', 'location', 'Northwest')
title('Aggregate Production, Consumption and Costs (logs)')
xlabel('Year')
xlim([beg_year, beg_year+Tmax])
saveas(gcf, 'figures/aggregates.png')

figure(3)
plot(periods, log(nrofpatents), '-k')
title('Number of patents produced by year')
xlim([beg_year, beg_year+Tmax])
xlabel('Year')
ylabel('Number of patents (logs)')
saveas(gcf, 'figures/num_pats.png')

figure(4)
plot(patent_shares(:,1), patent_shares(:,2), '-k', 'LineWidth', 1)
hold on
plot(patent_shares(:,1),patent_shares(:,3), '--k', 'LineWidth', 1)
plot(patent_shares(:,1),patent_shares(:,4), '-^k', 'LineWidth', 1, 'MarkerSize', 3) 
xlim([beg_year-1, beg_year+Tmax])
xlabel('Year')
ylabel('Patent Shares')
legend('New technologies', 'New combinations', 'Refinements')
title('Fraction of patents by nature of innovation')
hold off
saveas(gcf, 'figures/patents.png')

% Welfare gains with subsidy

welfareNT = squeeze(welfare_subs(:,:,1:length(subsidy_vals)));
welfareNC = squeeze(welfare_subs(:,:,length(subsidy_vals)+1:2*length(subsidy_vals)));
welfareBOTH = squeeze(welfare_subs(:,:,2*length(subsidy_vals)+1:end));

Delta_welfareNT = (welfareNT - repmat(welfare,1,length(subsidy_vals)))./repmat(welfare,1,length(subsidy_vals));
Delta_welfareNC = (welfareNC - repmat(welfare,1,length(subsidy_vals)))./repmat(welfare,1,length(subsidy_vals));
Delta_welfareBOTH = (welfareBOTH - repmat(welfare,1,length(subsidy_vals)))./repmat(welfare,1,length(subsidy_vals));

figure(5)
plot([0,subsidy_vals], [0,Delta_welfareNT(2,:)], '-k', [0,subsidy_vals], [0,Delta_welfareNC(2,:)], ...
    '--k', [0,subsidy_vals], [0,Delta_welfareBOTH(2,:)], '-^k')
title('Welfare Gains from Subsidy')
ylabel('Welfare gains (percentage terms, relative to no subsidy)')
xlabel('Subsidy Value')
legend('New Technologies', 'New Combinations', 'Both', 'location', 'Northwest')
saveas(gcf, 'figures/welfare_gains.png')




