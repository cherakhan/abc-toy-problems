% This script breaks down the differences between the Bayesian
% metropolis-hastings acceptance probability in textbooks and papers, and
% the form in takes when it is implemented in algorithms. In books it
% written as p(y|theta')*p(theta')/p(y|theta'')*p(theta''). However in
% algorithms in term is evaluated as it's negative natural log *2 and
% compressed into 1 exponential term. 

clear;clc

%% Likelihood and prior values

% likelihood value at timestep 1
lt1 = 0.5;
% prior value at timestep 1
pt1 = 0.05;

% likelihood value at timestep 2
lt2 = 0.4;
% prior value at timestep 2
pt2 = 0.07;

%% Compute acceptance probability equation form

alphaE = (lt2*pt2)/(lt1*pt1);

%% Compute 2*negative natual log values for l and p

% negative log likelihood at timestep 1
nlt1 = -2*log(lt1);

% negative log prior at timestep 1
npt1 = -2*log(pt1);

% negative log likelihood at timestep 2
nlt2 = -2*log(lt2);

% negative log prior at timestep 2
npt2 = -2*log(pt2);

%% Compute acceptance probability in algorithm form

alphaA = exp(-0.5*(nlt2-nlt1)-0.5*(npt2-npt1));

%% Output finding to the command window

if alphaE == alphaA
    fprintf('The acceptance probabilities are the same, alpha = %2.2f \n',alphaE)
end