function [] = plotchain(chain)
% Function to create plots of the mcmc chains generated during my thesis
%
% INPUT:
% chain = The markov chain for one parameter

iterations = (1:1:length(chain))';
plot(iterations,chain,'k')