function [nlogweight] = abcbananafun(par,data,tolerance)
% function to find the 2*negative log likelihood for the banana
% distribution ABC parameter estimation problem. This is solving the same
% examples as Marko Laine did except with a likelihood approximation.
% 
% INPUT:
% parameters = [muX,muY,var1,var2,correlation,bananityPar1,bananitypar2]
% data....
%
% OUTPUT:
% nloglikelihood = 2*log(p(y|x,theta))
%
% Author: Tom Connell
% Date: March 2018

%% Forward Simulate
CM = [par(3),par(5)*sqrt(par(3)*par(4));par(5)*sqrt(par(3)*par(4)),par(4)];
R = chol(CM);
sims = randn(1000,2)*R;
sims(:,1) = sims(:,1)+par(1);
sims(:,2) = sims(:,2)+par(2);
%bananaSims = bananafun(sims,[par(6),par(7)],1);
bananaSims = bananafun(sims,[1,1],1);
%% Compute summaries

obsSummaries = bananasummaries(data);
simSummaries = bananasummaries(bananaSims);

%% Comute error

distance = abs(simSummaries-obsSummaries);

%% Compute weighting function

% nlogweight = uniformkernel(distance,tolerance);
nlogweight = gaussiankernel(distance,tolerance);