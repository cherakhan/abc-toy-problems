function [error] = computeerror(par,data)
% Banana problem compute error

rho = par(5)*sqrt(par(3)*par(4));
cm = [par(3),rho;rho,par(4)];
r = chol(cm);
sims = randn(1000,2)*r;
sims(:,1) = sims(:,1)+par(1);
sims(:,2) = sims(:,2)+par(2);
sims = bananafun(sims,[1,1],1);
simsum = bananasummaries(sims);
datasum = bananasummaries(data);
error = abs(simsum-datasum);