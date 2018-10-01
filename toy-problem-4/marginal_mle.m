function [mle] = marginal_mle(chain,MIN,MAX)
% function to intake markov chain and give an accurate MLE estimate using
% the kde function (kenel density estimation)

[~,density,mesh] = kde(chain,2^7,MIN,MAX);

[~,i] = max(density);

mle = mesh(i);