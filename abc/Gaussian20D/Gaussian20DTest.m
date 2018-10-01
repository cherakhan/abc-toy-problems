% Script to test ABC-MCMC vs MCMC w/likelihood for a 2-D Gaussian problem
% trying to find mu

clear;clc

addpath('~/thesis/dram/dram')
addpath('~/thesis/dram/dram/utils')
addpath('~/thesis/dram/abc/abcutils')

for iMethod = 1:2
    
   if iMethod==1
       forward = @(par,data) (par'-data)'*inv(eye(20))*(par'-data);
   else
       tolerance = ones(20,1);
       forward = @(par,data) abcgaussian(par,data,tolerance);
   end
   
   % Sampler method
    adaptint = 1000;
    drscale = 5;
    
    % Number of iterations
    nsimu = 100000;
    
    % Chain starting location
    %start = ones(20,1)*10;
    start = fliplr(1:1:20);
    
    qcov = eye(20);
    
    %% define passed structs
    
    clear model data params options
    
    model.ssfun = forward;
    
    data = (1:1:20)';
    
    params.par0 = start;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    
    %% Call to DRAM
    [results,chain] = dramrun(model,data,params,options);
    
    %% Plot
    
    figure
    subplot(1,3,1)
    histogram(chain(:,1))
    subplot(1,3,2)
    histogram(chain(:,10))
    subplot(1,3,3)
    histogram(chain(:,20))
    
    
end