function return_vars = simulate_path_zeta(Tmax, T_subs, params, subsidy_vec, seed)

    rng(seed)

    etaH = params.etaH; 
    etaM = params.etaM; 
    etaL = params.etaL; 
    ttau = params.ttau; 
    llambda = params.llambda; 
    kkappa = params.kkappa; 
    xxi = params.xxi;   
    zeta = params.zeta; 
    ggamma = params.ggamma; 
    epsilon = params.epsilon; 
    rr = params.rr; 
    bbeta = 1/(1+rr); 
    PPi = ggamma/(1-bbeta); 
    nu = params.nu; 
    nrofinv = params.nrofinv; 
    nroffirms = round(nu*nrofinv); 
    g1 = params.g1; 
    g2 = params.g2; 

    year = 1836; % initial year
    M0 = 100; % initial number of technologies

    % Preallocate space for Mmat:
    blocksize = 50*M0;

    % [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
    Mmat = zeros(blocksize,5);
    Mmat(1:M0,1) = 1:M0;    % m ID
    Mmat(1:M0,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
    Mmat(1:M0,3) = zeros(M0,1);   % reuse/refinement so far
    Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
    Mmat(1:M0,5) = ones(M0,1); % technologies are born hot

    quality = ones(nroffirms,1); % initial quality of patents = 1
    nrofpatents = zeros(Tmax, 1);
    Growth = zeros(Tmax, 1);
    GDP = zeros(Tmax, 1);
    inventor_cost = zeros(Tmax,1);
    firm_cost = zeros(Tmax,1);
    policy_cost = zeros(Tmax,1);
    Consumption = zeros(Tmax,1);
    Discounted_utility = zeros(Tmax, 1);
    Discounted_utility_log = zeros(Tmax, 1);

    summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
    summat(1,:) = [year, 1, 0, 0]; %[year, new tech, new comb, reuse]
    seed = 1702;

    %%%%%%%%%%%%% DEFINE STATES
    %3 states: 1-new tech 2-recomb 3-reuse
    state = ones(1,M0);

    for t = 1:Tmax  % number of years of iteration. Here 200 years

        oldstate = state;
        year = year + 1;
        Mt = length(state);
        remove_line = [];

        % add new block of memory to Mmat if needed
        if( Mt + (blocksize/20) > length(Mmat) )         % less than 10% of blocksize free slots
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
        subsidy_cost = zeros(Mt,1);
        
        mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
        switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^( Mmat(1:Mt,2) - (year - Mmat(1:Mt,5))*zeta - 2)));   % (= 1 reuse, = 0 new comb)
        	% instead of 2^(n-2), I'm taking the expected value of n after t periods: n - t*zeta. t is how many periods the tech is hot

        if t < T_subs
            subsidy = [0,0];
        else
            subsidy = subsidy_vec;
        end

        for j=1:Mt

            % Draw new ideas and new costs
            % [mstar,aalpha] = modeldraws(j, ttau, llambda, kkappa); % not using seed
            y = rand(2,1); 
            mstar = ttau*log( (1 + y(1)*exp(1.0*j/ttau)) / (1-y(1)) );                     % Inverse of G
            aalpha = llambda*(- log(1 - y(2)) )^(1/kkappa) ;                           % Inverse of weibull CDF
            aalpha1 = aalpha*(1 - subsidy(1)); % new tech cost
            aalpha2 = aalpha*(1 - subsidy(2)); % new comb cost
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
                    elseif (mstar <= Mt && (aalpha2 < PPi*(etaM - etaL) ||  switchtocold(j) == 0  ) ) || ...
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
            patent_cost(j) = mmu(j)*(aalpha1*newtech + aalpha2*newcomb*(1-switchtocold(j)) ); % if product line if hot, new comb is free

            % cost of providing a subsidy for inventor's costs
            subsidy_cost(j) = mmu(j)*aalpha*( subsidy(1)*newtech + subsidy(2)*newcomb*(1-switchtocold(j)) ) ;
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
        policy_cost(t) = sum(subsidy_cost*mean(quality));
        
        Consumption(t) = GDP(t) - inventor_cost(t) - firm_cost(t) - policy_cost(t);
        if t > T_subs
            Discounted_utility(t) = bbeta^(t-T_subs-1) * Consumption(t)^(1-epsilon) / (1-epsilon) ;
            Discounted_utility_log(t) = bbeta^(t-T_subs-1) * log( Consumption(t) ) ;
        end

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
    Welfare_log = sum(Discounted_utility_log);

    return_vars = v2struct(Mmat, summat, GDP, Consumption, inventor_cost, firm_cost, policy_cost, nrofpatents, quality, Growth, Welfare, Welfare_log);

end









