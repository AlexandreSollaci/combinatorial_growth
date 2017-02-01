function [M1, M2, M3, M4] = moments(etaH, etaM, ttau, xxi)
    
%    etaH = .5; % NT step size
%    etaM = .1; % NC step size
    etaL = 0; % Refinement step size
    rr = 0.05; % interest rate
    RR = 1/(1+rr); % gross interest rate
    bbeta = 0.1; % labor share
    LL = 1; % production worker
    ppi = LL*(1-bbeta)*((bbeta/(1-bbeta))^bbeta);
    PPi = ppi/(1-RR);
%    ttau = 100; % shape parameter
    llambda = 1; % alpha has weibull distribution with scale parameter lambda and shape parameter kapa
    kkapa = .1;


    %%%%%%%%%%%%%%%%%%%%%%%%%%
%    xxi = 2;   %fraction of combos happening (if increased, peak shifts right)
    %low xxi increases the role of recomb e.g., xxi=0.1
    %high xxi lowers the role of recomb e.g., xxi=1000
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    year = 1836;
    M0 = 100; %initial number of technologies

    Mmat = zeros(M0,5);  % [ID, n (pool size), z (number of times reused), birth year, number of periods that a tech remains hot]
    Mmat(:,1) = 1:M0;    % m ID
    Mmat(:,2) = randi([1 3],M0,1);  % average tech pool in 1836 was 2.
    Mmat(:,3) = Mmat(:,1).*0;   % reuse/refinement so far
    Mmat(:,4) = Mmat(:,1).*0 + year;    % birth year
    Mmat(:,5) = ones(M0,1); % technologies are born hot

    Tmax = 50;
    %patenttype = [];
    patenttype = 2*diag(rand(M0,1))*[ones(M0,1), zeros(M0,2)];
    summat = zeros(Tmax + 1,5);
    summat(1,:) = [year, 0, sum(patenttype)/sum(sum(patenttype))];
    seed = 1702;

    nrofinv = 500; % number of inventors

    %%%%%%%%%%%%% DEFINE STATES
    %3 states: 1-new tech 2-recomb 3-reuse
    state = Mmat(:,1).*0 + 1;
    meann = zeros(Tmax,1);
    meanhotduration = zeros(Tmax,1);
    tic
    for t = 1:Tmax  % number of years of iteration. Here 200 years

        oldstate = state;
        year = year + 1;
        Mt = length(Mmat(:,1));
        mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line

        switchtocold = ( rand(Mt,1) < mmu.*xxi./(2.^(Mmat(:,2) - 2))) ;   % (= 1 reuse, = 0 new comb)

        for j=1:Mt

            [mstar,aalpha] = modeldraws(j, ttau, llambda, kkapa, seed); % assign the new draw
            newtech = 0;
            newcomb = 0;
            reuse = 0;

            if switchtocold(j) == 1 % whenever line becomes cold, restart count
                Mmat(j,5) = 0;
            else
                Mmat(j,5) = Mmat(j,5) + 1;
            end

            if oldstate(j)<=2  % used to be new technology or new combo

                if mstar > Mt && aalpha < switchtocold(j)*PPi*(etaH - etaL) + (1-switchtocold(j))*PPi*(etaH - etaM) %new tech happens
                    state(j) = 1;
    %                 %#####################
    %                 Mmat(j,2) = Mmat(j,2) + 1; %%%%%%%% when do new tech, bring in new tech to old pool??
    %                 %#####################
    %                 Mmat(j,3) = 0;
    %                 Mmat = [Mmat ; Mmat(j,:)]; 
    %                 Mmat(end,4) = year; 
    %                 Mmat(end,1) = length(Mmat(:,1));
                    Mmat = [Mmat ; [length(Mmat) + 1 , Mmat(j,2) + 1 , 0 , year, 1] ];
                    state = [state; 1];
                    newtech = 1;
                elseif (mstar <= Mt && (aalpha < PPi*(etaM - etaL) ||  switchtocold(j) == 0  )) || ...
                       (mstar > Mt && (1 - switchtocold(j))*(aalpha > PPi*(etaH - etaM)) > 0)    % new comb import
                    state(j) = 2;
                    Mmat(j,2) = Mmat(j,2) + 1;
                    Mmat(j,3) = 0;
                    newcomb = 1;
                else
                    state(j) = switchtocold(j)*3 + (1-switchtocold(j))*2;
                    Mmat(j,3) = switchtocold(j)*(Mmat(j,3) + 1) + (1-switchtocold(j))*0;
                    reuse = 1;
                end

            else %  oldstate(j) == 3   %used to be refinement
                if mstar > Mt && aalpha < PPi*(etaH - etaL) % new tech happens
                    state(j) = 1;
    %                 %##################### same as above. does a NT increase n?
    %                 Mmat(j,2) = Mmat(j,2)+1;
    %                 Mmat(j,3) = 0; %%%%%%%%%%%%% old tech line is still refinement, no?
    %                 Mmat = [Mmat;Mmat(j,:)]; 
    %                 Mmat(end,4) = year; 
    %                 Mmat(end,1) = length(Mmat(:,1));
                    Mmat = [Mmat ; [length(Mmat) + 1 , Mmat(j,2) + 1 , 0 , year, 1] ];
                    state = [state; 1];
                    newtech = 1;
                elseif mstar <= Mt && aalpha < PPi*(etaM - etaL) % new comb import
                    state(j) = 2;
                    Mmat(j,2) = Mmat(j,2)+1;
                    Mmat(j,3) = 0;
                    newcomb = 1;
                else
                    state(j) = 3;
                    Mmat(j,3) = Mmat(j,3)+1;
                    reuse = 1;
                end
            end
            patenttype = [patenttype; mmu(j)*newtech, mmu(j)*newcomb, mmu(j)*reuse]; % record patent types in each product line and year
        end

        nrofnewtech = sum(patenttype(:,1));
        nrofnewcomb = sum(patenttype(:,2));
        nrofreuse = sum(patenttype(:,3));

    % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % %    % if t>190
    % % % %         nn=Mmat(:,2)-(state<=2);
    % % % %         maxxx=max(nn)
    % % % %         newmat=zeros(length(Mmat(:,1)),maxxx);
    % % % %         for k=1:length(Mmat(:,1));
    % % % %             for l=1:nn(k);
    % % % %                 newmat(k,l)=nchoosek(nn(k),l);
    % % % %             end
    % % % %         end
    % % % %         meann=[];
    % % % %         for k=1:length(Mmat(:,1));
    % % % %             meann(k,1)=sum([1:1:maxxx].*newmat(k,:))/sum(newmat(k,:));
    % % % %         end
    % % % %         meann=meann+(state<=2);
    % % % %         
    % % % %         
    % % % %         summat(t,2)=mean(meann);
    % % % %     %end
    % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        meann(t) = mean(Mmat(Mmat(:,5)>0, 2));
        meanhotduration(t) = mean(Mmat(Mmat(:,5)>0,5));


        summat(t+1,1)= year;

        summat(t+1,3)= nrofnewtech/(nrofnewtech+nrofnewcomb+nrofreuse);
        summat(t+1,4)= nrofnewcomb/(nrofnewtech+nrofnewcomb+nrofreuse);
        summat(t+1,5)= nrofreuse/(nrofnewtech+nrofnewcomb+nrofreuse);

%        display(['Iteration number ' num2str(t)])
    end

    
    index1 = find(summat(:,1) == 1880);
    index2 = find(summat(:,1) == 1870);
    
    M1 = summat(index1, 5);
    M2 = summat(index1, 3);
    M3 = summat(index2, 5);
    M4 = summat(index2, 3);
    
end