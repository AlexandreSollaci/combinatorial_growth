clear
close all
clc

% turn on just in time compilation (makes for loops faster)
feature accel on 

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial_growth/simulations/simulate_model/')
rng(10)
etaH = 0.15; % NT step size
etaM = 0.1; % NC step size
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

year = 1836;
M0 = 100; %initial number of technologies

% Preallocate space for Mmat:
blocksize = 50*M0;

% [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
Mmat = zeros(blocksize,5);
Mmat(1:M0,1) = 1:M0;    % m ID
Mmat(1:M0,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
Mmat(1:M0,3) = zeros(M0,1);   % reuse/refinement so far
Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
Mmat(1:M0,5) = ones(M0,1); % technologies are born hot

Tmax = 180;

nu = 3.1; % chosen so that nrof firms makes avg growth rate is approximately 2%
nrofinv = exp(6); % number of inventors, chosen to match initial number of patents
nroffirms = round(nu*nrofinv); 
g1 = 0.066; % initial growth rate of patent numbers
g2 = 0.02; % final growth rate of patents

quality = ones(nroffirms,1); % initial quality of patents = 1
nrofpatents = zeros(Tmax, 1);
Growth = zeros(Tmax, 1);
GDP = zeros(Tmax, 1);
inventor_cost = zeros(Tmax,1);
firm_cost = zeros(Tmax,1);
Consumption = zeros(Tmax,1);
Discounted_utility = zeros(Tmax, 1);

summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
summat(1,:) = [year, 1, 0, 0]; %[year, new tech, new comb, reuse]

%%%%%%%%%%%%% DEFINE STATES
%3 states: 1-new tech 2-recomb 3-reuse
state = ones(1,M0);

for t = 1:Tmax  % number of years of iteration. Here 200 years

    oldstate = state;
    year = year + 1;
    Mt = length(state);
    remove_line = [];

    % add new block of memory to Mmat if needed
    if( Mt + (blocksize/10) > length(Mmat) )         % less than 10% of blocksize free slots
        Mmat(length(Mmat)+1:length(Mmat)+blocksize, :) = zeros(blocksize,5);       % add new BLOCK_SIZE slots
    end
    
    if year < 1896
    	nrofinv = (1 + g1)*nrofinv; % population grows at 6% rate
    	nroffirms = round(nu*nrofinv);
    else
    	nrofinv = (1 + g2)*nrofinv; % population grows at 2% rate
    	nroffirms = round(nu*nrofinv);
    end

    % Initialze variables
    growth = zeros(Mt,1);
    patent_type = zeros(Mt,3); % in each period, store the new patents that were created
    patent_cost = zeros(Mt,1); % keep track of inventors' costs
    
    mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
    switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^( Mmat(1:Mt,2) - (year - Mmat(1:Mt,5))*zeta - 2)));   % (= 1 reuse, = 0 new comb)
    	% instead of 2^(n-2), I'm taking the expected value of n after t periods: n - t*zeta. t is how many periods the tech is hot

    for j=1:Mt

        % Draw new ideas and new costs
        %[mstar,aalpha] = modeldraws(j, ttau, llambda, kkappa); % not using seed
        y = rand(2,1); 
        mstar = ttau*log( (1 + y(1)*exp(phi*j/ttau)) / (1-y(1)) );                     % Inverse of G
        aalpha = llambda*(- log(1 - y(2)) )^(1/kkappa) ;                               % Inverse of weibull CDF
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

        if rand(1) < zeta * Mmat(j, 2)
            Mmat(j,2) = max(Mmat(j,2) - 1, 0); % with prob n*zeta, tech pool looses one tech
        end

        if rand(1) < zeta
            remove_line = [remove_line, j]; % variable to store identity of dead lines
        else
            if oldstate(j)<=2  % used to be new technology or new combination

                % new technology
                if mstar > Mt && aalpha1 < switchtocold(j)*PPi*(etaH - etaL) + (1-switchtocold(j))*PPi*(etaH - etaM)
                    state(j) = 1;
                    aux = length(Mmat(Mmat(:,1)>0));
                    Mmat(aux + 1,:) = [aux + 1 , Mmat(j,2) + 1 , 0 , year, 1];
                    state = [state, 1];
                    newtech = 1; 

                % new combinations
                elseif (mstar <= Mt && (aalpha2 < PPi*(etaM - etaL) ||  switchtocold(j) == 0  )) || ...
                       (mstar > Mt && (1 - switchtocold(j))*(aalpha2 > PPi*(etaH - etaM)) > 0) 
                    state(j) = 2;
                    Mmat(j,2) = Mmat(j,2) + 1;
                    Mmat(j,3) = 0;
                    newcomb = 1;
                % reuse
                else
                    state(j) = switchtocold(j)*3 + (1-switchtocold(j))*2;
                    Mmat(j,3) = switchtocold(j)*(Mmat(j,3) + 1) + (1-switchtocold(j))*0;
                    reuse = 1;
                end

            else %  oldstate(j) == 3   % ==> used to be refinement

                % new technology
                if mstar > Mt && aalpha1 < PPi*(etaH - etaL) 
                    state(j) = 1;
                    aux = length(Mmat(Mmat(:,1)>0));
                    Mmat(aux + 1,:) = [aux + 1 , Mmat(j,2) + 1 , 0 , year, 1];
                    state = [state, 1];
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
        end
        
        % growth rate on product line j
        growth(j) = mmu(j)*( etaH*newtech + etaM*newcomb + etaL*reuse );

        % patent type in product line j, this year
        patent_type(j,:) = mmu(j)*[newtech, newcomb, reuse];

        % cost of producing patents this year in product line j
        patent_cost(j) = mmu(j)*aalpha*(newtech + newcomb);

    end

    % remove product lines that died
    Mmat(remove_line,:) = []; % with prob zeta, the product line dies
    state(remove_line) = [];

    % rate of growth of the economy in year t
    Growth(t) = sum(growth)/nroffirms;

    % aggregate variables for final good market clearing
    GDP(t) = sum(quality)/(1-ggamma) ;
    inventor_cost(t) = sum(patent_cost*mean(quality));
    firm_cost(t) = (1-ggamma)*sum(quality);
    
    Consumption(t) = GDP(t) - inventor_cost(t) - firm_cost(t);
    Discounted_utility(t) = bbeta^t * Consumption(t)^(1-epsilon) / (1-epsilon) ;

    % compute total number of patents by category
    nrofnewtech = sum(patent_type(:,1));
    nrofnewcomb = sum(patent_type(:,2));
    nrofreuse = sum(patent_type(:,3));
    
    % keep track of number of patents created each year
    nrofpatents(t) = nrofnewtech+nrofnewcomb+nrofreuse;
    
    % Populate summary matrix
    summat(t+1,1)= year;
    summat(t+1,2)= nrofnewtech/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new tech
    summat(t+1,3)= nrofnewcomb/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new combo
    summat(t+1,4)= nrofreuse/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of reuse

    % match patents with firms
    inventions = [etaH*ones(round(nrofnewtech),1); etaM*ones(round(nrofnewcomb),1); ...
        etaL*ones(round(nrofreuse),1)]; % vector with all new patents and their qualities
    inventions_comp = [inventions; zeros(nroffirms - length(inventions), 1)]; % complete vector with non-matched firms
    match = inventions_comp(randperm(length(inventions_comp))); % randomize inventions vector to get match to firms

    % compute the quality increase generated by patents this period
    quality_new = mean(quality)*match + [quality; zeros(length(match) - length(quality),1)]; % new quality for firms
        % note new firms start with quality = 0, otherwise we would be creating more patents out of nothing
    quality = quality_new; % set quality vector to new quality

    display(['Iteration number ' num2str(t)])
end

% get rid of any extra space allocated to MMat
Mmat = Mmat(Mmat(:,1)>0,:);

% Welfare 
Welfare = sum(Discounted_utility);

periods = 1836 + linspace(1,Tmax,Tmax);
figure(1)
plot(periods, Growth)
xlim([1836, year])
%ylim([0.015, 0.025])
title('Rate of Growth')
saveas(gcf, 'tex_files/figures/growth.png')

figure(2)
plot(periods, log(GDP), 'k', periods, log(Consumption), 'b', periods, log(firm_cost), 'r', periods, log(inventor_cost), 'g')
legend('GDP', 'Consumption', 'Firms` Cost', 'Inventors` Cost', 'location', 'Northwest')
title('Aggregate Production, Consumption and Costs (logs)')
xlabel('Year')
xlim([1836, year])
saveas(gcf, 'tex_files/figures/aggregates.png')

figure(3)
plot(periods, log(nrofpatents))
title('Number of patents produced by year')
xlim([1836, year])
xlabel('Year')
ylabel('Number of patents (logs)')
saveas(gcf, 'tex_files/figures/num_pats.png')

figure(4)
plot(summat(:,1),summat(:,2), 'r', summat(:,1),summat(:,3), 'g', summat(:,1),summat(:,4)) 
xlim([1836, year])
xlabel('Year')
ylabel('Patent Type Fractions')
legend('New technologies', 'New combinations', 'Refinements')
title('Fraction of patents by nature of innovation')
saveas(gcf, 'tex_files/figures/patents.png')

%%% Compute moments

[M, F] = compute_moments(summat);

display(F)

M1 = M(1);
M2 = M(2);
M3 = M(3);
M4 = M(4);
M5 = M(5);
M6 = M(6);
M7 = M(7);
M8 = M(8);

if Tmax < 110
    T = table;
    T.moment = {'model'; 'data'};
    T.new_tech1880 = [M1 ; 0.10];
    T.new_comb1880 = [M2 ; 0.35];
    T.reuse1880 = [M3 ; 0.55];
    T.new_tech1930 = [M4 ; 0.05];
    T.new_combo1930 = [M5 ; 0.5];
    T.reuse1930 = [M6 ; 0.4];

    T2 = table;
    T2.moment = {'model'; 'data'};
    T2.reuse_peak = [M7 ; 0.55];
    T2.peak_year = [M8 ; 34];

    display(T)
    display(T2)

else
    
    M9 = M(9);
    M10 = M(10);
    M11 = M(11);
    M12 = M(12);
    M13 = M(13);
    M14 = M(14);
    
    T = table;
    T.moment = {'model'; 'data'};
    T.new_tech1850 = [M1 ; 0.4];
    T.new_comb1850 = [M2 ; 0.25];
    T.reuse1850 = [M3 ; 0.35];
    T.new_tech1900 = [M4 ; 0.03];
    T.new_combo1900 = [M5 ; 0.45];
    T.reuse1900 = [M6 ; 0.52];

    T2 = table;
    T2.moment = {'model'; 'data'};
    T2.new_tech1950 = [M7 ; 0.02];
    T2.new_comb1950 = [M8 ; 0.75];
    T2.reuse1950 = [M9 ; 0.23];
    T2.new_tech2000 = [M10 ; 0.01];
    T2.new_comb2000 = [M11 ; 0.80];
    T2.reuse2000 = [M12 ; 0.19];

    T3 = table;
    T3.moment = {'model'; 'data'};
    T3.reuse_peak = [M13; 0.55];
    %T3.peak_year = [M14 ; 34];

    display(T)
    display(T2)
    display(T3)

    data = [0.4; 0.25; 0.35; 0.03; 0.45; 0.52; 0.02; 0.75; 0.33; 0.01; 0.80; 0.19; 0.55];% 34];
    Summary_stats = [data, M(1:end-1)'];

    rowLabels = {'new tech 1850'; 'new comb 1850'; 'reuse 1850'; 'new tech 1900'; 'new comb 1900'; 'reuse 1900'; ...
        'new tech 1950'; 'new comb 1950'; 'reuse 1950'; 'new tech 2000'; 'new comb 2000'; 'reuse 2000'; 'reuse peak'};% 'peak year'};
    columnLabels = {'Data', 'Model'};
    matrix2latex(Summary_stats, 'tex_files/moments_table.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels)
end














