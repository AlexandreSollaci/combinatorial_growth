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
ttau = 300; % shape parameter
llambda = 2.2; % alpha has weibull distribution with scale parameter lambda and shape parameter kappa
kkapa = 5;

ggamma = .5; % GUESS for gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%
xxi = 10;   %fraction of combos happening (if increased, peak shifts right)
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

Tmax = 100;
nrofinv = 500; % number of inventors

mmu0 = rand(M0,1)*2*nrofinv/M0; % initial measure of inventors in each product line
patenttype = [ mmu0, zeros(M0,2) ]; % keep track of the number of patents from each type: [new tech, new comb, reuse]
%patenttype = 2*diag(rand(M0,1))*[ones(M0,1), zeros(M0,2)];
summat = zeros(Tmax + 1,5); % matrix to keep track of fractions of each patent type
summat(1,:) = [year, 0, sum(patenttype)/sum(sum(patenttype))];
seed = 1702;

%%%%%%%%%%%%% DEFINE STATES
%3 states: 1-new tech 2-recomb 3-reuse
state = Mmat(:,1).*0 + 1;
meann = zeros(Tmax,1);
meanhotduration = zeros(Tmax,1);
Growth = zeros(Tmax, 1);

nroffirms = 660; % chosen so that avg growth rate is approximately 2%
quality = ones(nroffirms,1); % initial quality of patents

GDP = zeros(Tmax, 1);
inventor_cost = zeros(Tmax,1);
firm_cost = zeros(Tmax,1);

for t = 1:Tmax  % number of years of iteration. Here 200 years

    oldstate = state;
    year = year + 1;
    Mt = length(Mmat(:,1));
    mmu = rand(Mt,1)*2*nrofinv/Mt;  % number of inventors in product line
    
    growth = zeros(Mt,1);
    new_patent = zeros(Mt,3); % in each period, store the new patents that were created
    new_patent_cost = zeros(Mt,1); % keep track of inventors' costs

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
                %%#####################
                %Mmat(j,2) = Mmat(j,2) + 1; %%%%%%%% when do new tech, bring in new tech to old pool??
                %%#####################
                %Mmat(j,3) = 0;
                %Mmat = [Mmat ; Mmat(j,:)]; 
                %Mmat(end,4) = year; 
                %Mmat(end,1) = length(Mmat(:,1));
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
                %%##################### same as above. does a NT increase n?
                %Mmat(j,2) = Mmat(j,2)+1;
                %Mmat(j,3) = 0; %%%%%%%%%%%%% old tech line is still refinement, no?
                %Mmat = [Mmat;Mmat(j,:)]; 
                %Mmat(end,4) = year; 
                %Mmat(end,1) = length(Mmat(:,1));
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
        % record ALL patent types in each product line and year
        patenttype = [patenttype; mmu(j)*newtech, mmu(j)*newcomb, mmu(j)*reuse];
        
        % growth rate
        growth(j) = mmu(j)*( etaH*newtech + etaM*newcomb + etaL*reuse );

        % patent type in product line j, this year
        new_patent(j,:) = mmu(j)*[newtech, newcomb, reuse];

        % cost of producing patents this year
        new_patent_cost(j) = mmu(j)*aalpha;
    end

    % rate of growth of the economy in year t
    Growth(t) = sum(growth)/nroffirms; 

    % compute total number of patents by category
    pat_newtech = sum(new_patent(:,1));
    pat_newcomb = sum(new_patent(:,2));
    pat_reuse = sum(new_patent(:,3));
    
    % match patents with firms
    inventions = [etaH*ones(round(pat_newtech),1); etaM*ones(round(pat_newcomb),1); ...
        etaL*ones(round(pat_reuse),1)]; % vector with all new patents and their qualities
    inventions_aux = [inventions; zeros(nroffirms - length(inventions), 1)]; % complete vector with non-matched firms
    match = inventions_aux(randperm(length(inventions_aux))); % randomize inventions vector to get match to firms
    
    % compute the quality of products produced in current period
    quality_new = mean(quality)*match + quality; % new quality for firms
    quality = quality_new; % set quality vector to new quality

    % aggregate variables for goods market clear
    GDP(t) = sum(quality_new)/(1-ggamma) ;
    inventor_cost(t) = sum(new_patent_cost);
    firm_cost(t) = (1-ggamma)*sum(quality_new);
    

    nrofnewtech = sum(patenttype(:,1));
    nrofnewcomb = sum(patenttype(:,2));
    nrofreuse = sum(patenttype(:,3));

    summat(t+1,1)= year;

    summat(t+1,3)= nrofnewtech/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new tech
    summat(t+1,4)= nrofnewcomb/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of new combo
    summat(t+1,5)= nrofreuse/(nrofnewtech+nrofnewcomb+nrofreuse); % fraction of reuse

    display(['Iteration number ' num2str(t)])

    meann(t) = mean(Mmat(Mmat(:,5)>0, 2));
    meanhotduration(t) = mean(Mmat(Mmat(:,5)>0,5));
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % if t>190
    % % % %     nn=Mmat(:,2)-(state<=2);
    % % % %     maxxx=max(nn)
    % % % %     newmat=zeros(length(Mmat(:,1)),maxxx);
    % % % %     for k=1:length(Mmat(:,1));
    % % % %         for l=1:nn(k);
    % % % %             newmat(k,l)=nchoosek(nn(k),l);
    % % % %         end
    % % % %     end
    % % % %     meann=[];
    % % % % for k=1:length(Mmat(:,1));
    % % % %     meann(k,1)=sum([1:1:maxxx].*newmat(k,:))/sum(newmat(k,:));
    % % % % end
    % % % % meann=meann+(state<=2);
    % % % %         
    % % % %         
    % % % % summat(t,2)=mean(meann);
    % % % % %end
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

figure(1)
plot(1836:1935, Growth)
xlim([1836, 1935])
%ylim([0.015, 0.025])
title('Rate of Growth')

figure(2)
plot(summat(:,1),summat(:,3), 'r', summat(:,1),summat(:,4), 'g', summat(:,1),summat(:,5)) 
xlim([1836, year])
legend('New technologies', 'New combinations', 'Refinements')
title('Fraction of patents by nature of innovation')

%%% Compute moments

index1 = find(summat(:,1) == 1880);
index2 = find(summat(:,1) == 1930);

% M1 = summat(index1, 5); % = .55
M2 = summat(index1, 4);
M3 = summat(index1, 3);
% M4 = summat(index2, 5); % = /35
M5 = summat(index2, 4);
M6 = summat(index2, 3);
[M7 , M8] = max(summat(:,5)); % moments targetting the peak of reuse fraction and the year of the peak
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

% %% ######################## trying new set of moments ################
% [M1 , M2] = max(summat(:,5)); % moments targetting the peak of reuse fraction and the year of the peak
% [~, M3] = min( abs(summat(:,3) - summat(:,5)) ); % year where new tech and reuse cross
% [~, M4] = min( abs(summat(:,3) - summat(:,4)) ); % year where new tech and new comb cross
% M5 = summat(M3, 3); % value at which new tech and reuse cross
% M6 = summat(M4, 3);
% 
% F = (M1/0.6 - 1)^2 + (M2/34 - 1)^2 + (M3/24 - 1)^2 + (M4/29 - 1)^2 + (M5/0.37 - 1)^2 + (M6/.27 - 1)^2
% 
% T = table;
% T.moment = {'model'; 'data'};
% T.M1 = [M1 ; .60];
% T.M2 = [M2 ; 34];
% T.M3 = [M3 ; 24];
% T.M4 = [M4 ; 29];
% T.M5 = [M5 ; .37];
% T.M6 = [M6 ; .27];
    
% % figure(1)
% % subplot(2,2,1)
% % plot(summat(:,1),summat(:,3))
% % xlim([1836, year])
% % title('Fraction of New Technologies')
% % subplot(2,2,2)
% % plot(summat(:,1),summat(:,4))
% % xlim([1836, year])
% % title('Fraction of Novel Combinations')
% % subplot(2,2,3)
% % plot(summat(:,1),summat(:,5))
% % xlim([1836, year])
% % title('Fraction of Reuse/Refinement')
% % subplot(2,2,4)
% % plot(meann,meanhotduration)
% % title('Mean duration of hot line as function of mean n')
% % % plot(summat(:,1)+1836,summat(:,2))
% % % xlim([1836, year])
% % % title('Average Number of Technologies on a Patent')


% figure(1)
% plot(summat(:,1)+1800,summat(:,3),'LineWidth',1.75)
% title('Fraction of New Technologies','FontSize', 14)
% %xtitle('me')
% figure(2)
% plot(summat(:,1)+1800,summat(:,4),'LineWidth',1.75)
% title('Fraction of Novel Combinations','FontSize', 14)
% figure(3)
% plot(summat(:,1)+1800,summat(:,5),'LineWidth',1.75)
% title('Fraction of Reuse/Refinement','FontSize', 14)
% figure(4)
% plot(summat(:,1)+1800,summat(:,2),'LineWidth',1.75)
% title('Average Number of Technologies on a Patent','FontSize', 14)


