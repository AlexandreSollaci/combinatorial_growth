function  return_vars = simulate_path(Tmax, params, subsidy, quality, Mmat, summat_input, state)


    % return_vars = get_welfare(Tmax, params, subsidy, Mmat, summat, state)
    % where return_vars is a data struct with the following components:
    %   [Consumption, GDP, inventor_cost, firm_cost, policy_cost, Discounted_utility, Growth, Welfare, quality, Mmat, summat], state
    % Only Tmax and params are necessary to run; if the other arguments are not specified, the model starts from
    % period 0, that is 1836 (first year of patents).
    %
    % This function takes as input:
    %
    % Tmax: scalar, the number of periods over which the economy is to be interated on;
    %
    % subsidy: 2x1 vector, the subsidies to inventor's costs. This subsidy is computed as a percentage of
    % the cost, i.e. the cost will be multiplied by 1 - subsidy.
    % subsidy(1): subsidy for new technologies;
    % subsidy(2): subsidy for new combinations.
    %
    % params: the parameters for the simulation of the model.
    % IMPORTANT: params is a data struct!!

    etaH = params.etaH; 
    etaM = params.etaM;
    etaL = params.etaL;
    ttau = params.ttau;
    llambda = params.llambda;
    kkappa = params.kkappa;
    xxi = params.xxi;

    ggamma = params.ggamma;
    epsilon = params.epsilon;

    rr = params.rr;
	bbeta = 1/(1+rr); 
	PPi = ggamma/(1-bbeta); 

    nrofinv = params.nrofinv;
    nroffirms = params.nroffirms;

	M0 = 100; %initial number of technologies
	blocksize = 20*M0; % Preallocate space for Mmat

	if nargin == 2

		year = 1836; % initial year

		% [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
		Mmat = zeros(blocksize,5);
		Mmat(1:M0,1) = 1:M0;    % m ID
		Mmat(1:M0,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
		Mmat(1:M0,3) = zeros(M0,1);   % reuse/refinement so far
		Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
		Mmat(1:M0,5) = ones(M0,1); % technologies are born hot

		%3 states: 1-new tech 2-recomb 3-reuse
		state = ones(1,M0);
		quality = ones(nroffirms,1); % initial quality of patents = 1

		summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
		summat(1,:) = [year, 1, 0, 0]; %[year, new tech, new comb, reuse]

		subsidy = [0, 0];
	else

		year = summat_input(end, 1);
		summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
		summat(1,:) = summat_input(end, :); %[year, new tech, new comb, reuse]
		MMat = [Mmat ; zeros(blocksize,size(Mmat,2))];

	end


	Growth = zeros(Tmax, 1);
	GDP = zeros(Tmax, 1);
	inventor_cost = zeros(Tmax,1);
	firm_cost = zeros(Tmax,1);
	policy_cost = zeros(Tmax,1);
	Consumption = zeros(Tmax,1);
	Discounted_utility = zeros(Tmax, 1);

	for t = 1:Tmax  % number of years of iteration.

	    oldstate = state;
	    year = year + 1;
	    Mt = length(Mmat(Mmat(:,1)>0));
	    mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
	    
	    % add new block of memory to Mmat if needed
	    if( Mt + (blocksize/20) > length(Mmat) )         % less than 10% of blocksize free slots
	        Mmat(length(Mmat)+1:length(Mmat)+blocksize, :) = zeros(blocksize,5);       % add new BLOCK_SIZE slots
	    end

	    % Initialze variables
	    growth = zeros(Mt,1);
	    patent_type = zeros(Mt,3); % in each period, store the new patents that were created
	    patent_cost = zeros(Mt,1); % keep track of inventors' costs
	    subsidy_cost = zeros(Mt,1);
	    switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^(Mmat(1:Mt,2) - 2)) ) ;   % (= 1 reuse, = 0 new comb)

	    for j=1:Mt

	        % Draw new ideas and new costs
	        [mstar,aalpha] = modeldraws(j, ttau, llambda, kkappa); % not using seed
	        aalpha1 = aalpha * (1 - subsidy(1)); % new tech cost
	        aalpha2 = aalpha * (1 - subsidy(2)); % new comb cost
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
	        
	        % growth rate
	        growth(j) = mmu(j)*( etaH*newtech + etaM*newcomb + etaL*reuse );

	        % patent type in product line j, this year
	        patent_type(j,:) = mmu(j)*[newtech, newcomb, reuse];

	        % cost of producing patents this year
	        patent_cost(j) = mmu(j)*aalpha1*newtech + mmu(j)*aalpha2*newcomb*(switchtocold(j) == 1);

	        % subsidy cost
	        subsidy_cost(j) = mmu(j)*aalpha*( subsidy(1)*newtech + subsidy(2)*newcomb*(switchtocold(j)==1) );
	        	% note that new combinations are only costly when product line is cold

	    end

	    % rate of growth of the economy in year t
	    Growth(t) = sum(growth)/nroffirms;

	    % aggregate variables for final good market clearing
	    GDP(t) = sum(quality)/(1-ggamma) ;
	    inventor_cost(t) = sum(patent_cost*mean(quality));
	    firm_cost(t) = (1-ggamma)*sum(quality);
	    policy_cost(t) = sum(subsidy_cost*mean(quality));
	    Consumption(t) = GDP(t) - inventor_cost(t) - firm_cost(t) - policy_cost(t);
	    Discounted_utility(t) = bbeta^t * Consumption(t)^(1-epsilon) / (1-epsilon) ;

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

	% get rid of any extra space allocated to MMat
	Mmat = Mmat(Mmat(:,1)>0,:);

	% Welfare 
	Welfare = sum(Discounted_utility);

	return_vars = v2struct(Consumption, GDP, inventor_cost, firm_cost, policy_cost, Discounted_utility, ...
			Growth, Welfare, quality, Mmat, summat, state);

end