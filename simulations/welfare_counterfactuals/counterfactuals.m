clear
close all
clc

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/welfare_counterfactuals/code/')

etaH = 0.2; % NT step size
etaM = 0.05; % NC step size
etaL = 0; % Refinement step size
ttau = 800; % shape parameter for ideas distribution
llambda = 2; % alpha (cost) has weibull distribution with scale parameter lambda ...
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

Consumption1 = result1.Consumption;
GDP1 = result1.GDP;
inventor_cost1 = result1.inventor_cost;
firm_cost1 = result1.firm_cost;
policy_cost1 = result1.policy_cost;
Discounted_utility1 = result1.Discounted_utility;
Growth1 = result1.Growth;
Welfare1 = result1.Welfare;
quality1 = result1.quality;
summat1 = result1.summat;
state1 = result1.state;

Mmat1 = result1.Mmat;
numtech1 = zeros(T_past, 1);

for t = 1:T_past
    index = find(Mmat1(:,4) == 1836 + t, 1, 'last');
    numtech1(t) = Mmat1(index, 1);
end

% define parameters for counterfactuals
T_forward = 50;
subsidy = [0, 0; 0.25, 0]; % sibsidy for new tech, subsidy for new combo

% initialize matrices to store results
Consumption2 = zeros(T_forward, size(subsidy,1));
GDP2 = zeros(T_forward, size(subsidy,1));
inventor_cost2 = zeros(T_forward, size(subsidy,1));
firm_cost2 = zeros(T_forward, size(subsidy,1));
policy_cost2 = zeros(T_forward, size(subsidy,1));
Discounted_utility2 = zeros(T_forward, size(subsidy,1));
Growth2 = zeros(T_forward, size(subsidy,1));
numtech2 = zeros(T_forward, size(subsidy,1));
Welfare2 = zeros(1, size(subsidy,1));
summat2 = zeros(T_forward + 1, 4,  size(subsidy,1));

% run simulation for counterfactual
for i = 1:size(subsidy,1)
    result2 = simulate_path(T_forward, params, subsidy(i,:), quality1, Mmat1, summat1, state1);

    Consumption2(:,i) = result2.Consumption;
    GDP2(:,i) = result2.GDP;
    inventor_cost2(:,i) = result2.inventor_cost;
    firm_cost2(:,i) = result2.firm_cost; 
    policy_cost2(:,i) = result2.policy_cost;
    Discounted_utility2(:,i) = result2.Discounted_utility;
    Growth2(:,i) = result2.Growth; 
    Welfare2(i) = result2.Welfare; 
    summat2(:,:,i) = result2.summat; % [year, new tech, new combo, reuse]

    Mmat2 = result2.Mmat;

    for t = 1:T_forward
        index = find(Mmat2(:,4) == Mmat1(end,4) + t, 1, 'last');
        numtech2(t, i) = length(Mmat2(1:index, 1));
    end

end

Consumption = [repmat(Consumption1, 1,  size(subsidy,1)); Consumption2];
GDP = [repmat(GDP1, 1,  size(subsidy,1)); GDP2];
inventor_cost = [repmat(inventor_cost1, 1,  size(subsidy,1)); inventor_cost2];
firm_cost = [repmat(firm_cost1, 1,  size(subsidy,1)); firm_cost2];
policy_cost = [repmat(policy_cost1, 1,  size(subsidy,1)); policy_cost2];
Discounted_utility = [repmat(Discounted_utility1, 1,  size(subsidy,1)); Discounted_utility2];
Growth = [repmat(Growth1, 1,  size(subsidy,1)); Growth2];
Welfare = [repmat(Welfare1, 1,  size(subsidy,1)); Welfare2];
summat = [repmat(summat1(:,:), 1, 1, size(subsidy, 1)) ; summat2(2:end,:,:)];
numtech = [repmat(numtech1, 1,  size(subsidy,1)); numtech2];

Tmax = T_past + T_forward;
periods = 1836 + linspace(1,Tmax,Tmax);
periods2 = 1836 + T_past + linspace(1, T_forward, T_forward);

aggregates = figure(1);
subplot(3,2,1)
plot(periods, log(Consumption(:,1)), 'r', periods, log(Consumption(:,2)), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Consumption')
xlabel('Year'); ylabel('Consumption')

subplot(3,2,2)
plot(periods, log(GDP(:,1)), 'r', periods, log(GDP(:,2)), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of GDP')
xlabel('Year'); ylabel('GDP')

subplot(3,2,3)
plot(periods, log(firm_cost(:,1)), 'r', periods, log(firm_cost(:,2)), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Firm costs')
xlabel('Year'); ylabel('Firms` Costs')

subplot(3,2,4)
plot(periods, log(inventor_cost(:,1)), 'r', periods, log(inventor_cost(:,2)), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Inventor costs')
xlabel('Year') ;ylabel('Inventors` Costs')

subplot(3,2,5)
plot(periods, log(policy_cost(:,1)), 'r', periods, log(policy_cost(:,2)), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Policy costs')
xlabel('Year') ;ylabel('Subsidy Costs')

subplot(3,2,6)
plot(periods, Growth(:,1), 'r', periods, Growth(:,2), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northeast')
title('Growth rate of the economy')
xlabel('Year'); ylabel('Growth rate')
xlim([1836, 1836+Tmax])
saveas(gcf, 'aggregates.png')

aggregates2016 = figure(2);
subplot(3,2,1)
plot(periods2, log(Consumption2(:,1)), 'r', periods2, log(Consumption2(:,2)), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Consumption, starting in 2016')
xlabel('Year'); ylabel('Consumption')

subplot(3,2,2)
plot(periods2, log(GDP2(:,1)), 'r', periods2, log(GDP2(:,2)), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of GDP, starting in 2016')
xlabel('Year'); ylabel('GDP')

subplot(3,2,3)
plot(periods2, log(firm_cost2(:,1)), 'r', periods2, log(firm_cost2(:,2)), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Firm costs, starting in 2016')
xlabel('Year'); ylabel('Firms` Costs')

subplot(3,2,4)
plot(periods2, log(inventor_cost2(:,1)), 'r', periods2, log(inventor_cost2(:,2)), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Inventor costs, starting in 2016')
xlabel('Year') ;ylabel('Inventors` Costs')

subplot(3,2,5)
plot(periods2, log(policy_cost2(:,1)), 'r', periods2, log(policy_cost2(:,2)), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('log of Policy costs, starting in 2016')
xlabel('Year') ;ylabel('Subsidy Costs')

subplot(3,2,6)
plot(periods2, Growth2(:,1), 'r', periods2, Growth2(:,2), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Growth rate of the economy, starting in 2016')
xlim([1836+T_past, 1836+Tmax])
saveas(gcf, 'aggregates2016.png')

aggregates_diff = figure(3);
subplot(3,1,1)
plot(periods2, log(Consumption2(:,2)) - log(Consumption2(:,1)))
xlim([1836+T_past, 1836+Tmax])
title('log diff of Consumption, with and without subsidy')
xlabel('Year'); ylabel('log difference of Consumption')

subplot(3,1,2)
plot(periods2, log(GDP2(:,2)) - log(GDP2(:,1)))
xlim([1836+T_past, 1836+Tmax])
title('log diff of GDP, with and without subsidy')
xlabel('Year'); ylabel('log difference of GDP')

subplot(3,1,3)
plot(periods2, log(firm_cost2(:,2)) - log(firm_cost2(:,1)))
xlim([1836+T_past, 1836+Tmax])
title('log diff of Firms costs, with and without subsidy')
xlabel('Year'); ylabel('log difference of Firms costs')
saveas(gcf, 'aggregates_diff.png')

newtechcreation = figure(4);
subplot(2,1,1)
plot(periods(2:end), diff(numtech(:,1)), 'r', periods(2:end), diff(numtech(:,2)), '--b')
legend('No subsidy', 'Subsdidy')
title('Number of technologies created')
xlabel('Year'); ylabel('Number of new techs')
xlim([1836, 1836+Tmax])

subplot(2,1,2)
plot(periods2(2:end), diff(numtech2(:,1)), 'r', periods2(2:end), diff(numtech2(:,2)), '--b')
legend('No subsidy', 'Subsdidy')
title('Number of technologies created, starting in 2016')
xlabel('Year'); ylabel('Number of new techs')
xlim([1836+T_past, 1836+Tmax])
saveas(gcf, 'tech_creation.png')

patents = figure(5);
subplot(2,2,1)
plot(summat(:,1,1), summat(:,2,1), 'r', summat(:,1,2), summat(:,2,2), '--b')
legend('No subsidy', 'Subsdidy')
title('New technology evolution')
xlabel('Year'); ylabel('Fraction of new techs')
xlim([1836, 1836+Tmax])

subplot(2,2,2)
plot(summat(:,1,1), summat(:,3,1), 'r', summat(:,1,2), summat(:,3,2), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Southeast')
title('New combination evolution')
xlabel('Year'); ylabel('Fraction of new combos')
xlim([1836, 1836+Tmax])

subplot(2,2,3)
plot(summat(:,1,1), summat(:,4,1), 'r', summat(:,1,2), summat(:,4,2), '--b')
legend('No subsidy', 'Subsdidy')
title('Reuse/refinement evolution')
xlabel('Year'); ylabel('Fraction of reuse')
xlim([1836, 1836+Tmax])

subplot(2,2,4)
plot(periods, numtech(:,1), 'r', periods, numtech(:,2), '--b')
xlim([1836, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Number of existing technologies')
xlabel('Year') ;ylabel('Number of product lines')
saveas(gcf, 'patents.png')

patents2016 = figure(6);
subplot(2,2,1)
plot(summat2(:,1,1), summat2(:,2,1), 'r', summat2(:,1,2), summat2(:,2,2), '--b')
legend('No subsidy', 'Subsdidy')
title('New technology evolution, starting in 2016')
xlabel('Year'); ylabel('Fraction of new techs')
xlim([1836+T_past, 1836+Tmax])

subplot(2,2,2)
plot(summat2(:,1,1), summat2(:,3,1), 'r', summat2(:,1,2), summat2(:,3,2), '--b')
legend('No subsidy', 'Subsdidy', 'location', 'Southeast')
title('New combination evolution, starting in 2016')
xlabel('Year'); ylabel('Fraction of new combos')
xlim([1836+T_past, 1836+Tmax])

subplot(2,2,3)
plot(summat2(:,1,1), summat2(:,4,1), 'r', summat2(:,1,2), summat2(:,4,2), '--b')
legend('No subsidy', 'Subsdidy')
title('Reuse/refinement evolution, starting in 2016')
xlabel('Year'); ylabel('Fraction of reuse')
xlim([1836+T_past, 1836+Tmax])

subplot(2,2,4)
plot(periods2, numtech2(:,1), 'r', periods2, numtech2(:,2), '--b')
xlim([1836+T_past, 1836+Tmax])
legend('No subsidy', 'Subsdidy', 'location', 'Northwest')
title('Number of existing technologies, starting in 2016')
xlabel('Year') ;ylabel('Number of product lines')
saveas(gcf, 'patents2016.png')

T = table;
T.Welfare_no_sub = Welfare2(1);
T.Welfare_sub = Welfare2(2);
T.Welfare_increase = sign(Welfare2(1))*( Welfare2(2) - Welfare2(1) ) / Welfare2(1);
T.Numtech_no_sub = numtech(end,1);
T.Numtech_sub = numtech(end,2);
T.Numtech_increase = (numtech(end,2) - numtech(end,1)) / numtech(end,1);

T2 = table;
T2.Average_growth_no_sub = mean(Growth2(:,1));
T2.Average_growth_sub = mean(Growth2(:,2));
T2.Average_growth_increase = ( mean(Growth2(:,2)) - mean(Growth2(:,1)) )/ mean(Growth2(:,1));

display(T)
display(T2)


Summary_stats = [Welfare2 , sign(Welfare2(1))*( Welfare2(2) - Welfare2(1) ) / Welfare2(1) ; ...
    numtech(end,:) , (numtech(end,2) - numtech(end,1)) / numtech(end,1); ...
    mean(Growth2) , ( mean(Growth2(:,2)) - mean(Growth2(:,1)) )/ mean(Growth2(:,1)) ];

rowLabels = {'Welfare'; '\# technologies'; 'Average GDP growth'};
columnLabels = {'No Subsidy', 'Subsidy', '\% Change'};
matrix2latex(Summary_stats, 'counterfactual_table.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels)

