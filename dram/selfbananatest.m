% Test of using DRAM, same as bananatest.m
% Author: Tom Connell
% Date: March 2018

clear;clc

% '/home/tom/thesis/dram/utils'
addpath([pwd,filesep,'utils']);

%% Define input variables

method = 'DRAM';

% Number of iterations
nsimu = 10000;

% Scale for second proposal
drscale = 2; 

% how often to adapt
adaptint = 10; 

% chain starting location + data fed into forward
mu = [0,0];

% 'target covariance', fed into forward as data
cmat = [1, 0.9; 0.9, 1]; 
imat = inv(cmat);

% 'Bananity' of the target, fed into forward as data
bpar = [1,1]; 

% The forward problem. Making this function call returns -2log(likelihood)
% i.e the sum-of-squares/lest-squares error
bananass = @(x,d) bananafun(x-d.mu,d.bpar,1)*d.imat*bananafun(x-d.mu,d.bpar,1)';

%[xmin,ymin;xmax,ymax]
bounds = [-10,-10;10,10];


%% Define passed structs

clear model data params options

model.ssfun = bananass;

data = struct('mu',mu,'imat',imat,'bpar',bpar);

params.par0 = mu;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = eye(2)*5;


%% Call to dramrun

[results,chain] = dramrun(model,data,params,options);


%% Plotting

plotbivariate(chain(:,1),chain(:,2),'','')
figure
subplot(2,1,1)
plotchain(chain(:,1))
subplot(2,1,2)
plotchain(chain(:,2))