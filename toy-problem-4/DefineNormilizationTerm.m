% This is a script to define the term which will normilize the variance of
% each distribution of abc-residual summary statistics to 1. This is in
% effect taking into account the way in which summaries varies across
% different scales. It can also act as a Monte Carlo simulator for the
% posterior with a few extra steps at the end

addpath('../dram')
addpath('../abc/abcutils')


%% Monte Carlo simulate parameters values

% kSims Monte Carlo simulations in the range:
% Gaussian Pars muX -5:5, muY -5:5, varX 0:10, varY 0:10, correlation 0:1
% All Pars muX -5:5, muY -5:5, varX 0:10, varY 0:10, correlation 0:1,
% B1 0:5, B2 0:5

% kSims = 5000000;
kSims = 10000;

%% Gaussian Parameters

muX = (rand(kSims,1)*10)-5;
muY = (rand(kSims,1)*10)-5;
varX = rand(kSims,1)*10;
varY = rand(kSims,1)*10;
correlation = rand(kSims,1);

B1 = ones(kSims,1);
B2 = ones(kSims,1);


%% All Parameters

% muX = (rand(kSims,1)*10)-5;
% muY = (rand(kSims,1)*10)-5;
% varX = rand(kSims,1)*10;
% varY = rand(kSims,1)*10;
% correlation = (rand(kSims,1)*2)-1;
% B1 = rand(kSims,1)*5;
% B2 = rand(kSims,1)*5;


%% Forward simulate data + record residual

% container for kSims x 8 summaries and residuals
% summaries: MC simulation vector of statistics
% residuals: absolute value distance between obsSummaries and simSummaries
summaries = zeros(kSims,8);
residuals = zeros(kSims,8);

obs = csvread('bananaData.dat');
obsSums = bananasummaries(obs);

for ii = 1:kSims
    
    % Set up simulation variables
    mu = [muX(ii);muY(ii)];
    covar = correlation(ii)*sqrt(varX(ii)*varY(ii));
    Sigma = [varX(ii),covar;covar,varY(ii)];
    R = chol(Sigma);
    bananity = [B1(ii),B2(ii)];
    
    % Stochastically simulate
    sims = randn(100,2)*R;
    sims(:,1) = sims(:,1)+mu(1);
    sims(:,2) = sims(:,2)+mu(2);                        
    sims = bananafun(sims,bananity,1);
    
    % Compute summaries and abs-distance
    summaries(ii,:) = bananasummaries(sims);
    residuals(ii,:) = abs(bananasummaries(sims)-obsSums);
    
    
end


%% Define and output the Normilization term (vector)

normTerm = [ std(residuals(:,1)), std(residuals(:,2)), std(residuals(:,3)),...
    std(residuals(:,4)), std(residuals(:,5)), std(residuals(:,6)),...
    std(residuals(:,7)), std(residuals(:,8))];

% Gaussian pars
csvwrite('normTermGaussPars.csv',normTerm);

% All pars
%csvwrite('normTermAllPars.csv',normTerm);


%% Bonus rejection sampler for posterior

toRejectOrNotToReject = 0;

if toRejectOrNotToReject
    
    % ??? Not sure if needed
    % Single matrix of MC simulations
%     simParamsBanana = [muX,muY,varX,varY,correlation];
    simParamsBanana = [muX,muY,varX,varY,correlation,B1,B2];

    % Normalize distribution of summary statistics
    % A.k.a Normalise residual variance to 1
    residuals(:,1) = residuals(:,1)./normTerm(1);
    residuals(:,2) = residuals(:,2)./normTerm(2);
    residuals(:,3) = residuals(:,3)./normTerm(3);
    residuals(:,4) = residuals(:,4)./normTerm(4);
    residuals(:,5) = residuals(:,1)./normTerm(5);
    residuals(:,6) = residuals(:,2)./normTerm(6);
    residuals(:,7) = residuals(:,3)./normTerm(7);
    residuals(:,8) = residuals(:,4)./normTerm(8);

    % tolerance vector
    tolerance = [0.005/10,0.025/10,0.025/10,0.025/10,0.025/10,0.025/10,0.025/10,0.025/10];

%     posterior = zeros(kSims,5);
    posterior = zeros(kSims,7);
    for ii = 1:kSims
        % acceptance probability
        sumOfSquares = gaussiankernel(residuals(ii,:),tolerance);
        alpha = exp(-0.5*sumOfSquares);
        if rand <= alpha
%             posterior(ii,:) = simParamsBanana(ii,:);
            posterior(ii,:) = simParamsBanana(ii,:);
        end
    end

    % Individual containers for posterior samples
    postMuX = posterior(:,1);
    postMuY = posterior(:,2);
    postVarX = posterior(:,3);
    postVarY = posterior(:,4);
    postCorr = posterior(:,5);
    postB1 = posterior(:,6);
    postB2 = posterior(:,7);

    % Remove 0 values
    postMuX = postMuX(postMuX~=0);
    postMuY = postMuY(postMuY~=0);
    postVarX = postVarX(postVarX~=0);
    postVarY = postVarY(postVarY~=0);
    postCorr = postCorr(postCorr~=0);
    postB1 = postB1(postB1~=0);
    postB2 = postB2(postB2~=0);

    % Plot marginal posterior samples (Gaussian pars)

    figure(1)

    subplot(4,2,1)
    histogram(postMuX)
    xlabel('\mu_X')

    subplot(4,2,2)
    histogram(postMuY)
    xlabel('\mu_Y')

    subplot(4,2,3)
    histogram(postVarX)
    xlabel('\sigma^2_X')

    subplot(4,2,4)
    histogram(postVarY)
    xlabel('\sigma^2_Y')

    subplot(4,2,5)
    histogram(postCorr)
    xlabel('\rho')
    
    subplot(4,2,6)
    histogram(postB1)
    xlabel('B1')
    
    subplot(4,2,7)
    histogram(postB2)
    xlabel('B2')
    

    % Realisations of joint-distribution (posterior)
    for ii = 1:length(postMuX)
        figure(2)
        scatter(obs(:,1),obs(:,2),'k.')
        hold on
        bananaplot(0,0,1,1,0.9,1,1,'k-')
        bananaplot(postMuX(ii),postMuY(ii),postVarX(ii),postVarY(ii),postCorr(ii),1,1,'r-')
        xlim([-10,10])
        ylim([-5,15])
        saveas(gcf,['out',filesep,'post',num2str(ii),'.png'])
        clf
    end

    posterior = [postMuX,postMuY,postVarX,postVarY,postCorr,postB1,postB2];

end