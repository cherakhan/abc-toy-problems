function [range] = normterm(model,data)
% This is a script to define the term which will normilize the variance of
% each distribution of abc-residual summary statistics to 1. This is in
% effect taking into account the way in which summaries varys across
% different scales. This is a function version of the tried and tested
% script DefineNormilizationTerm.m

addpath('../dram')
addpath('../abc/abcutils')


%% Compute observed summaries from data

obs_sums = bananasummaries(data);


%% Monte Carlo simulate parameters values

% kSims Monte Carlo simulations in the range:
% Gaussian Pars muX -5:5, muY -5:5, varX 0:10, varY 0:10, correlation 0:1
% All Pars muX -5:5, muY -5:5, varX 0:10, varY 0:10, correlation 0:1

kSims = 10000;


%% Gaussian Parameters

if strcmp('Gaussian',model)

    muX = (rand(kSims,1)*10)-5;
    muY = (rand(kSims,1)*10)-5;
    varX = rand(kSims,1)*10;
    varY = rand(kSims,1)*10;
    correlation = rand(kSims,1);
    b1 = ones(kSims,1);
    b2 = ones(kSims,1);
    
end


%% All Parameters

if strcmp('fullbanana',model)

    muX = (rand(kSims,1)*10)-5;
    muY = (rand(kSims,1)*10)-5;
    varX = rand(kSims,1)*10;
    varY = rand(kSims,1)*10;
    correlation = rand(kSims,1);
    b1 = rand(kSims,1)*5;
    b2 = rand(kSims,1)*5;
    
end


%% Forward simulate data + record residual

% Container for kSims x 8 summaries and residuals
% summaries: MC simulation vector of statistics
% residuals: absolute value distance between obs_summaries and sim_summaries
summaries = zeros(kSims,length(obs_sums));
residuals = zeros(kSims,length(obs_sums));

for ii = 1:kSims
    
    % Set up simulation variables
    mu = [muX(ii);muY(ii)];
    covar = correlation(ii)*sqrt(varX(ii)*varY(ii));
    Sigma = [varX(ii),covar;covar,varY(ii)];
    R = chol(Sigma);
    bananity = [b1(ii),b2(ii)];
    
    % Stochastically simulate
    sims = randn(length(data),2)*R;
    sims(:,1) = sims(:,1)+mu(1);
    sims(:,2) = sims(:,2)+mu(2);                        
    sims = bananafun(sims,bananity,1);
    
    % Compute summaries and abs-distance
    summaries(ii,:) = bananasummaries(sims);
    residuals(ii,:) = abs(bananasummaries(sims)-obs_sums);
    
end


%% Define and output the Normilization term (vector)

norm_term = [std(residuals(:,1)), std(residuals(:,2)), std(residuals(:,3)),...
    std(residuals(:,4)), std(residuals(:,5)), std(residuals(:,6)),...
    std(residuals(:,7)), std(residuals(:,8))];

range = norm_term;