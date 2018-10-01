function [nlogweight] = gaussiankernel(d,t)
% function to implement the gaussian weighting function
% for p(y|y*,theta). Returns the negative log for direct use in dram.
%
% INPUT: 
% d = vector of distances between the set of summary statistics
%       normally always absolute norm
% t = vector of tolerances, 1 for each summary statistics

dt = d./t;
nlogweight = sum(dt.^2);