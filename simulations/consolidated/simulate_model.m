function outputs = simulate_model(params, nrofperiods, beg_year, beg_inv, summary_matrix, patent_matrix,...
 tech_state, prod_quality, subsidy, seed)

    % turn on just in time compilation (makes for loops faster)
    feature accel on 

    rng(seed)

    %% Load parameters
    etaH = params.etaH;
    etaM = params.etaM;
    etaL = params.etaL;
    tau = params.tau;
    phi = params.phi;
    lambda = params.lambda;
    kappa = params.kappa;
    xi = params.xi;
    zeta = params.zeta;
    gamma = params.gamma;
    epsilon = params.epsilon;
    beta = params.beta;
    Pi = params.Pi;
    nu = params.nu;
    g1 = params.g1;
    g2 = params.g2;
    Tdur = params.Tdur;

    %% Load initial variables

    Tmax = nrofperiods;
    year = beg_year;
    nrofinv = beg_inv;
    nroffirms = round(nu*nrofinv); 

    Mmat = summary_matrix;
    patmat = patent_matrix;
    state = tech_state;
    quality = prod_quality;

    if subsidy == [0 0]
        subsidyNT = zeros(Tmax, 1);
        subsidyNC = zeros(Tmax, 1);
    else
        subsidyNT = [subsidy(1)*ones(Tdur, 1); zeros(Tmax - Tdur, 1)];
        subsidyNC = [subsidy(2)*ones(Tdur, 1); zeros(Tmax - Tdur, 1)];
    end

    blocksize = 5000; % preallocate space for matrices

    nrofpatents = zeros(Tmax, 1);
    growth = zeros(Tmax, 1);
    GDP = zeros(Tmax, 1);
    inventor_cost = zeros(Tmax,1);
    firm_cost = zeros(Tmax,1);
    policy_cost = zeros(Tmax,1);
    consumption = zeros(Tmax,1);
    discounted_utility = zeros(Tmax,1);
    avg_quality = zeros(Tmax,1);
    techline_size = zeros(Tmax,1);

    %% Start loops

    for t = 1:Tmax  % number of years of iteration.

        oldstate = state(state>0);
        year = year + 1;
        Mt = length(oldstate);
        remove_line = [];

        % add new block of memory to Mmat if needed
        if( Mt + (blocksize/10) > length(Mmat) )         % less than 10% of blocksize free slots
            Mmat(length(Mmat)+1:length(Mmat)+blocksize, :) = zeros(blocksize,5);       % add BLOCKSIZE new slots
            state(length(state)+1:length(state)+blocksize) = zeros(1,blocksize);
        end
        
        if year < 1896
        	nrofinv = (1 + g1)*nrofinv; % population grows at 6% rate
        	nroffirms = round(nu*nrofinv);
        else
        	nrofinv = (1 + g2)*nrofinv; % population grows at 2% rate
        	nroffirms = round(nu*nrofinv);
        end

        % Initialze variables
        innov_size = zeros(Mt,1);
        patent_type = zeros(Mt,3); % in each period, store the new patents that were created
        patent_cost = zeros(Mt,1); % keep track of inventors' costs
        subsidy_cost = zeros(Mt,1); % keep track of subsidy costs
        
        mu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
        % (= 1 line turns cold, = 0 line remains hot)
        switchtocold = ( rand(Mt,1) < mu.*xi./(2.^( Mmat(1:Mt,2) - (year - Mmat(1:Mt,5))*zeta - 2)));   
        	% instead of 2^(n-2), I'm taking the expected value of n after t periods: n - t*zeta. t is how many periods the tech is hot

        for j=1:Mt

            % Draw new ideas and new costs
            y = rand(2,1); 
            mstar = tau*log( (1 + y(1)*exp(phi*j/tau)) / (1-y(1)) );                     % Inverse of G
            alpha0 = lambda*(- log(1 - y(2)) )^(1/kappa) ;                               % Inverse of weibull CDF
            alphaNT = alpha0*(1-subsidyNT(t)); % new tech cost
            alphaNC = alpha0*(1-subsidyNC(t)); % new comb cost
            % cost of reuse normalized to zero

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

                % initalize patent type counters
                newtech = 0;
                newcomb = 0;
                reuse = 0;

                if oldstate(j)<=2  % used to be new technology or new combination

                    % new technology
                    if mstar > Mt && alphaNT <= switchtocold(j)*Pi*(etaH - etaL) + (1-switchtocold(j))*Pi*(etaH - etaM)
                        state(j) = 1;
                        aux = length(Mmat(Mmat(:,1)>0));
                        Mmat(aux+1,:) = [aux+1, Mmat(j,2)+1, 0, year, 1];
                        state(aux+1) = 1;
                        newtech = 1; 

                    % new combinations
                    elseif (mstar <= Mt && (alphaNC <= Pi*(etaM - etaL) ||  switchtocold(j) == 0  )) || ...
                       (mstar > Mt && switchtocold(j)==0 && alphaNT > Pi*(etaH - etaM) ) 
                        state(j) = 2;
                        Mmat(j,2) = Mmat(j,2) + 1;
                        Mmat(j,3) = 0;
                        newcomb = 1;

                    % reuse
                    else
                        state(j) = 3;
                        Mmat(j,3) = Mmat(j,3) + 1;
                        reuse = 1;
                    end

                else %  oldstate(j) == 3 ==> used to be refinement

                    % new technology
                    if mstar > Mt && alphaNT <= Pi*(etaH - etaL) 
                        state(j) = 1;
                        aux = length(Mmat(Mmat(:,1)>0));
                        Mmat(aux+1, :) = [aux + 1, Mmat(j,2)+1, 0, year, 1];
                        state(aux+1) = 1;
                        newtech = 1;

                    % new combination   
                    elseif mstar <= Mt && alphaNC <= Pi*(etaM - etaL) 
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
            
            
                % innovation size of product line j in year t
                innov_size(j) = mu(j)*( etaH*newtech + etaM*newcomb + etaL*reuse);

                % patent type in product line j in year t
                patent_type(j,:) = mu(j)*[newtech, newcomb, reuse];

                hotline = (oldstate(j) <= 2) * (switchtocold(j) == 0); % = 1 if product line is hot

                % cost of producing patents this year in product line j in year t
                patent_cost(j) = mu(j)*(alphaNT*newtech + alphaNC*newcomb*(1-hotline)); % if product line is hot, new comb is free

                % cost of providing a subsidy for inventor's costs in year t
                subsidy_cost(j) = mu(j)*alpha0*( subsidyNT(t)*newtech + subsidyNC(t)*newcomb*(1-hotline)) ;

            end

        end

        % remove product lines that died
        if length(remove_line) < length(Mmat(Mmat(:,1)>0))
            Mmat(remove_line,:) = [];
            state(remove_line) = [];
        else
            Mmat(remove_line(1:end-1),:) = [];
            state(remove_line(1:end-1)) = [];
        end

        % rate of growth of the economy in year t
        growth(t) = sum(innov_size)/nroffirms;

        % aggregate variables for final good market clearing
        GDP(t) = sum(quality)/(1-gamma) ;
        firm_cost(t) = (1-gamma)*sum(quality);
        inventor_cost(t) = sum(patent_cost*mean(quality));
        policy_cost(t) = sum(subsidy_cost*mean(quality));
        
        consumption(t) = GDP(t) - inventor_cost(t) - firm_cost(t) - policy_cost(t);
        %discounted_utility(t) = beta^t * (consumption(t)^(1-epsilon) - 1)/ (1-epsilon) ;
        discounted_utility(t) = beta^t * log(consumption(t)) ;

        % compute total number of patents by category
        nrofnewtech = sum(patent_type(:,1));
        nrofnewcomb = sum(patent_type(:,2));
        nrofreuse = sum(patent_type(:,3));
        
        % keep track of number of patents created each year
        nrofpatents(t) = nrofnewtech+nrofnewcomb+nrofreuse;
        
        % Populate patent matrix
        if nrofpatents(t) > 0
            patmat(t,1)= year;
            patmat(t,2)= nrofnewtech/nrofpatents(t);    % fraction of new tech
            patmat(t,3)= nrofnewcomb/nrofpatents(t);    % fraction of new combo
            patmat(t,4)= nrofreuse/nrofpatents(t);      % fraction of reuse
        elseif nrofpatents(t) == 0 && t > 1 
            patmat(t,:) = patmat(t-1,:);                % If all product lines die, no patents are created
        else
            patmat(t,:) = [year, 1, 0, 0];
        end
        
        % record average quality and technology line size
        avg_quality(t) = mean(quality);
        techline_size(t) = length(state(state>0));

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

    % get rid of any extra space allocated to Mmat
    Mmat = Mmat(Mmat(:,1)>0,:);
    state = state(state>0);

    % Welfare 
    welfare = sum(discounted_utility);

   outputs = v2struct(Mmat, state, patmat, GDP, consumption, inventor_cost, firm_cost, policy_cost, ...
    nrofpatents, quality, growth, welfare, nrofinv, avg_quality, techline_size);

end





