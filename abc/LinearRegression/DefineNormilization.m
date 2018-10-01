% Define the normilization term

k = 10000;

% Define parameter range (Use the same as will be used by the sampler)
slope = rand(k,1)*20;
intercept = (rand(k,1)*40)-20;

% Read in data
obsData = csvread('data.csv');
obsSummaries = linregsummaries(obsData);

% Containers
summaries = zeros(k,2);
residuals = zeros(k,2);

for ii = 1:k
    
    % forward simulate
    x = [(-50:1:50)',ones(101,1)];
    beta = [slope(ii);intercept(ii)];
    y = x*beta;
    y = y+(randn*5); % plus noise (same as observed data creation)
    
    % compute summaries
    simSummaries = linregsummaries(y);
    summaries(ii,:) = simSummaries;
    
    % compute distance
    distance = abs(simSummaries-obsSummaries);
    residuals(ii,:) = distance;
end

normTerm = [std(residuals(:,1)),std(residuals(:,2))];
csvwrite('normilizationTerm.csv',normTerm)