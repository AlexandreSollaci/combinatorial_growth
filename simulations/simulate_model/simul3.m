clear all
close all
clc

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/')

etaH = 0.15; % NT step size
etaM = 0.03; % NC step size
etaL = 0; % Refinement step size
rr = 0.05; % interest rate
RR = 1/(1+rr); % gross interest rate
bbeta = 0.1; % labor share
LL = 1; % production worker
ppi = LL*(1-bbeta)*((bbeta/(1-bbeta))^bbeta);
PPi = ppi/(1-RR);

ttau = 300; % shape parameter for ideas distribution
llambda = 2.2; % alpha (cost) has weibull distribution with scale parameter lambda ...
kkapa = 5; % and shape parameter kappa

ggamma = .8; % GUESS for gamma
epsilon = 2; % GUESS for epsilon

%%%%%%%%%%%%%%%%%%%%%%%%%%
xxi = 10;   %fraction of combos happening (if increased, peak shifts right)
%low xxi increases the role of recomb e.g., xxi=0.1
%high xxi lowers the role of recomb e.g., xxi=1000
%%%%%%%%%%%%%%%%%%%%%%%%%%

year = 1836;
M0 = 100; %initial number of technologies

% [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
Mmat = zeros(M0,5);
Mmat(:,1) = 1:M0;    % m ID
Mmat(:,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
Mmat(:,3) = Mmat(:,1).*0;   % reuse/refinement so far
Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
Mmat(:,5) = ones(M0,1); % technologies are born hot

Tmax = 100;
nrofinv = 500; % number of inventors
nroffirms = 660; % chosen so that avg growth rate is approximately 2%
quality = ones(nroffirms,1); % initial quality of patents = 1

Growth = zeros(Tmax, 1);
GDP = zeros(Tmax, 1);
inventor_cost = zeros(Tmax,1);
firm_cost = zeros(Tmax,1);
Consumption = zeros(Tmax,1);
Discounted_utility = zeros(Tmax, 1);

summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
summat(1,:) = [year, 1, 0, 0]; %[year, new tech, new comb, reuse]
seed = 1702;

%%%%%%%%%%%%% DEFINE STATES
%3 states: 1-new tech 2-recomb 3-reuse
state = Mmat(:,1).*0 + 1;

for t = 1:Tmax  % number of years of iteration. Here 200 years

    oldstate = state;
    year = year + 1;
    Mt = length(Mmat(:,1));
    mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
    
    growth = zeros(Mt,1);
    patent_type = zeros(Mt,3); % in each period, store the new patents that were created
    patent_cost = zeros(Mt,1); % keep track of inventors' costs

    switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^(Mmat(:,2) - 2))) ;   % (= 1 reuse, = 0 new comb)

    for j=1:Mt

        % Draw new ideas and new costs
        [mstar,aalpha] = modeldraws(j, ttau, llambda, kkapa, seed);
        aalpha1 = aalpha; % new tech cost
        aalpha2 = aalpha; % new comb cost
        % cost of reuse normalized to zero

        % initalize patent type counters
        newtech = 0;
        newcomb = 0;
        reuse = 0;

        if switchtocold(j) == 1 % whenever line becomes cold, restart count
            Mmat(j,5) = 0;
        else
            Mmat(j,5) = Mmat(j,5) + 1;
        end

        if oldstate(j)<=2  % used to be new technology or new combination

            % new technology
            if mstar > Mt && aalpha1 < switchtocold(j)*PPi*(etaH - etaL) + (1-switchtocold(j))*PPi*(etaH - etaM)
                state(j) = 1;
                Mmat = [Mmat ; [length(Mmat) + 1 , Mmat(j,2) + 1 , 0 , year, 1] ];
                state = [state; 1];
                newtech = 1; 

            % new combinations
            elseif (mstar <= Mt && (aalpha2 < PPi*(etaM - etaL) ||  switchtocold(j) == 0  )) || ...
                   (mstar > Mt && (1 - switchtocold(j))*(aalpha2 > PPi*(etaH - etaM)) > 0) 
                state(j) = 2;
                Mmat(j,2) = Mmat(j,2) + 1;
                Mmat(j,3) = 0;
                newcomb = 1;
            else
                state(j) = switchtocold(j)*3 + (1-switchtocold(j))*2;
                Mmat(j,3) = switchtocold(j)*(Mmat(j,3) + 1) + (1-switchtocold(j))*0;
                reuse = 1;
            end

        else %  oldstate(j) == 3   % ==> used to be refinement

            % new technology
            if mstar > Mt && aalpha1 < PPi*(etaH - etaL) 
                state(j) = 1;
                Mmat = [Mmat ; [length(Mmat) + 1 , Mmat(j,2) + 1 , 0 , year, 1] ];
                state = [state; 1];
                newtech = 1;

            % new combination   
            elseif mstar <= Mt && aalpha2 < PPi*(etaM - etaL) 
                state(j) = 2;
                Mmat(j,2) = Mmat(j,2)+1;
                Mmat(j,3) = 0;
                newcomb = 1;

            % reuse
            else
                state(j) = 3;
                Mmat(j,3) = Mmat(j,3)+1;
                reuse = 1;
            end
        end
        
        % growth rate
        growth(j) = mmu(j)*( etaH*newtech + etaM*newcomb + etaL*reuse );

        % patent type in product line j, this year
        patent_type(j,:) = mmu(j)*[newtech, newcomb, reuse];

        % cost of producing patents this year
        patent_cost(j) = mmu(j)*aalpha;
    end

    % rate of growth of the economy in year t
    Growth(t) = sum(growth)/nroffirms;

    % aggregate variables for final good market clearing
    GDP(t) = sum(quality)/(1-ggamma) ;
    inventor_cost(t) = sum(patent_cost*mean(quality));
    firm_cost(t) = (1-ggamma)*sum(quality);
    Consumption(t) = GDP(t) - inventor_cost(t) - firm_cost(t);
    Discounted_utility(t) = RR^t * Consumption(t)^(1-epsilon) / (1-epsilon) ;

    % compute total number of patents by category
    nrofnewtech = sum(patent_type(:,1));
    nrofnewcomb = sum(patent_type(:,2));
    nrofreuse = sum(patent_type(:,3));
    
    % match patents with firms
    inventions = [etaH*ones(round(nrofnewtech),1); etaM*ones(round(nrofnewcomb),1); ...
        etaL*ones(round(nrofreuse),1)]; % vector with all new patents and their qualities
    inventions_comp = [inventions; zeros(nroffirms - length(inventions), 1)]; % complete vector with non-matched firms
    match = inventions_comp(randperm(length(inventions_comp))); % randomize inventions vector to get match to firms

    % compute the quality increase generated by patents this period
    quality_new = mean(quality)*match + quality; % new quality for firms
    quality = quality_new; % set quality vector to new quality

    % Populate summary matrix
    summat(t+1,1)= year;
    summat(t+1,2)= nrofnewtech/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new tech
    summat(t+1,3)= nrofnewcomb/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new combo
    summat(t+1,4)= nrofreuse/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of reuse

    display(['Iteration number ' num2str(t)])
end

% Welfare 
Welfare = sum(Discounted_utility);

figure(1)
periods = 1836 + linspace(1,100,100);
plot(periods, Growth)
xlim([1836, year])
%ylim([0.015, 0.025])
title('Rate of Growth')

figure(2)
plot(periods, GDP, 'k', periods, Consumption, 'b', periods, firm_cost, 'r', periods, inventor_cost, 'g')
legend('GDP', 'Consumption', 'Firms` Cost', 'Inventors` Cost')
title('Aggregate Production, Consumption and Costs')
xlabel('Year')
xlim([1836, year])

figure(3)
plot(periods, Discounted_utility)
title('Discounted Utility')
xlim([1836, year])
xlabel('Year')
ylabel('Utility of Representative Consumer')

figure(4)
plot(summat(:,1),summat(:,2), 'r', summat(:,1),summat(:,3), 'g', summat(:,1),summat(:,4)) 
xlim([1836, year])
xlabel('Year')
ylabel('Patent Fractions')
legend('New technologies', 'New combinations', 'Refinements')
title('Fraction of patents by nature of innovation')


%%% Compute moments

index1 = find(summat(:,1) == 1880);
index2 = find(summat(:,1) == 1930);

% M1 = summat(index1, 5); % = .55
M2 = summat(index1, 3);
M3 = summat(index1, 2);
% M4 = summat(index2, 5); % = /35
M5 = summat(index2, 3);
M6 = summat(index2, 2);
[M7 , M8] = max(summat(:,4)); % moments targetting the peak of reuse fraction and the year of the peak
% peak is aroung 60% and on year 1870 (1870 - 1836 = 34). 

F = (M2 - 0.30)^2 + (M3 - .1)^2 + (M5 - 0.60)^2 + (M6 - .03)^2 + (M7 - 0.6)^2 

T = table;
T.moment = {'model'; 'data'};
%T.M1 = [M1 ; .55];
T.M2 = [M2 ; .30];
T.M3 = [M3 ; .10];
%T.M4 = [M4 ; .35];
T.M5 = [M5 ; .60];
T.M6 = [M6 ; .03];
T.M7 = [M7 ; .60];
T.M8 = [M8 ; 34]
