function [nloglikelihood] = abclikelihoodgauss(par,data,tolerance,kernel)
% Function to find the 2*negative log likelihood for the bivariate Gaussian
% parameter estimation problem for my thesis. 
%
% Summary statistics: mean x, mean y, variance x, variance y, covariance
% 
% INPUT:
% parameters = row vector of 5 parameters
%   (1. MuX,2. MuY,3. VarX,4.VarY, 5. Correlation)
% data = the observed data for the problem
%
% OUTPUT
% nloglikelihood = 2*log(p(y|x,theta))
%
% Author: Tom Connell
% Date: March 2018


%% Forward simulate

% define covariance
CV = par(5)*sqrt(par(3)*par(4));

% assemble covariance matrix
CM = [par(3),CV;CV,par(4)];

% cholsky factorization
Rsim = chol(CM);

% simulate Gaussian random samples
sims = (randn(100,2))*Rsim;

% distribute according to mu value
xSims = sims(:,1)+par(1);
ySims = sims(:,2)+par(2);

%% Compute summaries

% for the simulations
simsMuX = mean(xSims);
simsMuY = mean(ySims);
simsCov = cov(xSims,ySims);

summaries = [simsMuX;simsMuY;simsCov(1,1);simsCov(2,2);sims(1,2)];
distance = abs(summaries-data);

%% Uniform kernel

if kernel =='u'
    nloglikelihood = uniformkernel(distance,tolerance);
end

%% Gaussian kernel

if kernel =='g'
    nloglikelihood = gaussiankernel(distance,tolerance);
end

