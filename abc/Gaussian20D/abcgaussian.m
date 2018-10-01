function [nlogweight] = abcgaussian(par,data,tolerance)

%% forward simulate

nSims = 20;
sims = randn(nSims,20) + meshgrid(par);

%% compute summaries

obsSummaries = data;
simSummaries = mean(sims)';
distance = abs(obsSummaries-simSummaries);

%% compute weight

nlogweight = gaussiankernel(distance,tolerance);