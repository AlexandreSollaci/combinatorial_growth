function F = MMobj(guess, Tmax)

    etaH = guess.etaH; % NT step size
    etaM = guess.etaM; % NC step size
    etaL = 0; % Refinement step size
    ttau = guess.tau; % shape parameter for ideas distribution
    llambda = guess.lambda; % alpha (cost) has weibull distribution with scale parameter lambda ...
    kkappa = guess.kappa; % and shape parameter kappa
    xxi = guess.xxi;   % 1/xi is the fraction of feasible combinations

    nrofinv = 500; % number of inventors

    ggamma = .4; % match labor share
    epsilon = 2; % intertemporal substitution

    rr = 0.05; % interest rate
    bbeta = 1/(1+rr); % discount factor
    LL = 1; % production worker
    PPi = ggamma/(1-bbeta); % constant part of the price of innovation (proposition 1)

    year = 1836; % initial year
    M0 = 100; %initial number of technologies

    % Preallocate space for Mmat:
    blocksize = 20*M0;

    % [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
    Mmat = zeros(blocksize,5);
    Mmat(1:M0,1) = 1:M0;    % m ID
    Mmat(1:M0,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
    Mmat(1:M0,3) = zeros(M0,1);   % reuse/refinement so far
    Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
    Mmat(1:M0,5) = ones(M0,1); % technologies are born hot

    summat = zeros(Tmax + 1,4); % matrix to keep track of fractions of each patent type
    summat(1,:) = [year, 1, 0, 0]; %[year, new tech, new comb, reuse]
    seed = 1702;

    %%%%%%%%%%%%% DEFINE STATES
    %3 states: 1-new tech 2-recomb 3-reuse
    state = ones(1,M0);

    for t = 1:Tmax  % number of years of iteration. Here 200 years

        oldstate = state;
        year = year + 1;
        Mt = length(Mmat(Mmat(:,1)>0));
        mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
        
        % add new block of memory to Mmat if needed
        if( Mt + (blocksize/10) > length(Mmat) )         % less than 10% of blocksize free slots
            Mmat(length(Mmat)+1:length(Mmat)+blocksize, :) = zeros(blocksize,5);       % add new BLOCK_SIZE slots
        end

        % Initialze variables
        growth = zeros(Mt,1);
        patent_type = zeros(Mt,3); % in each period, store the new patents that were created
        patent_cost = zeros(Mt,1); % keep track of inventors' costs

        switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^(Mmat(1:Mt,2) - 2)) ) ;   % (= 1 reuse, = 0 new comb)

        for j=1:Mt

            % Draw new ideas and new costs
            [mstar,aalpha] = modeldraws(j, ttau, llambda, kkappa); % not using seed
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
            patent_cost(j) = mmu(j)*aalpha;

        end

        % rate of growth of the economy in year t multiplied by number of firms
        JGrowth(t) = sum(growth);

        % calibrate number of firms so that average growth equals 2%
        nroffirms = mean(JGrowth) / 0.02;

        % compute total number of patents by category
        nrofnewtech = sum(patent_type(:,1));
        nrofnewcomb = sum(patent_type(:,2));
        nrofreuse = sum(patent_type(:,3));

        % Populate summary matrix
        summat(t+1,1)= year;
        summat(t+1,2)= nrofnewtech/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new tech
        summat(t+1,3)= nrofnewcomb/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new combo
        summat(t+1,4)= nrofreuse/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of reuse

        %display(['Iteration number ' num2str(t)])
    end

    % get rid of any extra space allocated to MMat
    Mmat = Mmat(Mmat(:,1)>0,:);


    %%% Compute moments %%%

    index1 = find(summat(:,1) == 1880);
    index2 = find(summat(:,1) == 1930);

    M1 = summat(index1, 2);
    M2 = summat(index1, 3);
    M3 = summat(index1, 4); % = .55

    M4 = summat(index2, 2);
    M5 = summat(index2, 3);
    M6 = summat(index2, 4); % = /35
    [M7 , M8] = max(summat(:,4)); % moments targetting the peak of reuse fraction and the year of the peak
    % peak is aroung 60% and on year 1870 (1870 - 1836 = 34). 

    F = (M1/0.1 - 1)^2  + (M2/0.35 - 1)^2 + (M3/0.55 - 1)^2 + (M4/0.03 - 1)^2 + (M5/0.5 - 1)^2 + ...
    (M6/0.4 - 1)^2 + (M7/0.55 - 1)^2 + (M8/34 - 1)^2 ;

end