function [nlogweight] = uniformkernel(d,t)
% function to implement the uniform weighting function
% for p(y|y*,theta). Returns the negative log for direct use in dram.
%
% INPUT: 
% d = vector of distances between the set of summary statistics
%       normally always absolute norm
% t = vector of tolerances, 1 for each summary statistics

dt = d./t;
if dt <= 1
    likelihood = 1;
else
    likelihood = 0;
end
nlogweight = -2*log(likelihood);