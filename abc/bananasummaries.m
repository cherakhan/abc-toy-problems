function [summaries] = bananasummaries(data)
% function which intakes bananified data and returns summary statistics
% about the data

% All parameters
polySummaries = polyfit(data(:,1),data(:,2),2);
stdx = std(data(:,1));
stdy = std(data(:,2));
meanx = mean(data(:,1));
meany = mean(data(:,2));
skewnessy = skewness(data(:,2));

summaries = [polySummaries,stdx,stdy,meanx,meany,skewnessy];

% Mu only
% meanx = mean(data(:,1));
% meany = mean(data(:,2));
% summaries = [meanx,meany];
 
% Gaussian par only
% stdx = std(data(:,1));
% stdy = std(data(:,2));
% meanx = mean(data(:,1));
% meany = mean(data(:,2));
% summaries = [stdx,stdy,meanx,meany];