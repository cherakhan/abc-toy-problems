function [summaries] = linregsummaries(y)

% As of Vrugt
summaries = [mean(y);std(y)];