function [nlogweight] = abclikelihood(par,tolerance,data)

%% Forward simulate

x = [(-50:1:50)',ones(101,1)];
beta = [par(1);par(2)];
epsilon = randn(length(x),1)*par(3);

y = (x*beta) + epsilon;

%% Compute summaries

simSummaries = linregsummaries(y);
obsSummaries = linregsummaries(data);

%% Compute distance

distance = abs(simSummaries-obsSummaries);

%% Compute weighting kernel

nlogweight = uniformkernel(distance,tolerance);