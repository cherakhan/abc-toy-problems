addpath('~/thesis/dram/dram')

% kSims Monte Carlo simulations in the range:
% muX -5:5
% muY -5:5
% varX 0:10
% varY 0:10
% correlation -1:1
% B1 0:5
% B5 0:5

kSims = 10000;

muX = (rand(kSims,1)*10)-5;
muY = (rand(kSims,1)*10)-5;
varX = rand(kSims,1)*10;
varY = rand(kSims,1)*10;
correlation = (rand(kSims,1)*2)-1;
b1 = rand(kSims,1)*5;
b2 = rand(kSims,1)*5;

% container for kSims x 8 summaries and residuals
% summaries: MC simulation vector of statistics
% residuals: absolute value distance with each normalized to var = 1
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
    bananity = [b1(ii),b2(ii)];
    
    % Stochastically simulate
    sims = randn(100,2)*R;
    sims(:,1) = sims(:,1)+mu(1);
    sims(:,2) = sims(:,2)+mu(2);                        
    sims = bananafun(sims,bananity,1);
    
    % Compute summaries and abs-distance
    summaries(ii,:) = bananasummaries(sims);
    residuals(ii,:) = abs(bananasummaries(sims)-obsSums);
end

% Single matrix of MC simulations
simParamsBanana = [muX,muY,varX,varY,correlation,b1,b2];

% Normalise residual variance to 1
tuning = [std(residuals(:,1));std(residuals(:,2));std(residuals(:,3));std(residuals(:,4));std(residuals(:,5));std(residuals(:,6));std(residuals(:,7));std(residuals(:,8))];
csvwrite('tuning.csv',tuning)
residuals(:,1) = residuals(:,1)./std(residuals(:,1));
residuals(:,2) = residuals(:,2)./std(residuals(:,2));
residuals(:,3) = residuals(:,3)./std(residuals(:,3));
residuals(:,4) = residuals(:,4)./std(residuals(:,4));
residuals(:,5) = residuals(:,5)./std(residuals(:,5));
residuals(:,6) = residuals(:,6)./std(residuals(:,6));
residuals(:,7) = residuals(:,7)./std(residuals(:,7));
residuals(:,8) = residuals(:,8)./std(residuals(:,8));

% Implement standard rejection sampler (uniform kernel)
residuals = (sum(residuals'))';

tolerance = 0.5;

indicator = residuals<=tolerance;
% Output number of posterior samples
sum(indicator)
simParamsBanana = simParamsBanana.*indicator;

% Matrix to hold the posterior samples
postMuX = simParamsBanana(:,1);
postMuY = simParamsBanana(:,2);
postVarX = simParamsBanana(:,3);
postVarY = simParamsBanana(:,4);
postCorr = simParamsBanana(:,5);
postB1 = simParamsBanana(:,6);
postB2 = simParamsBanana(:,7);
posterior = [postMuX(postMuX~=0),postMuY(postMuY~=0),postVarX(postVarX~=0),postVarY(postVarY~=0),postCorr(postCorr~=0),postB1(postB1~=0),postB2(postB2~=0)];

% Plotting marginal posterior for each unknown
figure(1)
histogram(posterior(:,1));
figure(2)
histogram(posterior(:,2));
figure(3)
histogram(posterior(:,3));
figure(4)
histogram(posterior(:,4));
figure(5)
histogram(posterior(:,5));
figure(6)
histogram(posterior(:,6));
figure(7)
histogram(posterior(:,7));

% Loop to print figure for each joint posterior sample
% data = csvread('bananaData.dat');
%
% for ii = 1:sum(indicator)
%     figure(1)
%     scatter(data(:,1),data(:,2),'k.')
%     hold on
%     bananaplot(0,0,1,1,0.9,1,1,'k-')
%     bananaplot(sol(ii,1),sol(ii,2),sol(ii,3),sol(ii,4),sol(ii,5),sol(ii,6),sol(ii,7),'r-')
%     xlim([-10,10])
%     ylim([-5,15])
%     saveas(gcf,['out',filesep,'post',num2str(ii),'.png'])
%     clf
% end
