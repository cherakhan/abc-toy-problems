function [alpha] = acceptanceRule(a,b)
% Function to be called by the dram code to replace the problematic line
% where exp(-0.5*(Inf-Inf)) is Nan and not zero. 

if isnan(b)
    alpha = 0;
else
    alpha = min(a,b);
end