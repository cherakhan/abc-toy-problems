function [nloglikelihood] = analyticallikelihood(par,data)
% function to take in the model parameters m, b (for a straight line,
% y= mx+ b) and the standard deviation of the noise term sigma,
% to a normal dist. N(0,sigma^2). Then fit these parameters given the data
% with traditional likelihood machinery. The function returns 
% -2log(p(y|theta))
% m = par(1), b = par(2), sigma (sd) = par(3)

%% Forward simulate 

N = 101;

x = [(-50:1:50)',ones(101,1)];
beta = [par(1);par(2)];
y = x*beta;

%% Compute negative log likelihood

ss = sum((data-y).^2);

%% Normalize

loglikelihood = -N/2*log(2*pi) -N/2*log(par(3)^2) -(1/(2*par(3)^2))*ss;

%% Re-negate

nloglikelihood = -2*loglikelihood;