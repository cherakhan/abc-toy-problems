% Create noisy data from a linear model
% 101 samples from -50-50

x = [(-50:1:50)',ones(101,1)];
% [slope;intercept]
beta = [5;10];
% clean data
y = x*beta;
% noise 
noise = randn(length(y),1)*5;
% add noise to data
y = y+noise;

% output
csvwrite('data.csv',y)