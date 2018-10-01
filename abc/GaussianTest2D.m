% Analytical Gaussian test 2D

clear;clc

addpath([pwd,filesep,'utils']);

methods = {'mh','am','dr','dram'};

%% Run all methods
for iMethods = 1:4
    
    method = methods{iMethods};
    
    switch method
        case 'mh'
            adaptint = 0;
            drscale = 0;
            plotTitle = 'M-H';
        case 'am'
            adaptint = 100;
            drscale = 0;
            plotTitle = 'AM';
        case 'dr'
            adaptint = 0;
            drscale = 2;
            plotTitle = 'DR';
        case 'dram'
            adaptint = 100;
            drscale = 2;
            plotTitle = 'DRAM';
    end
    
    %% Define input variables
    
    % number of simulations
    nsimu = 1000;
    
    % starting location
    start = [5,5];
    
    % unbounded
    
    % returns -2log(likelihood)
    forward = @(par,data) (par-data.mu)*data.ICM*(par-data.mu)';
    
    % proposal distribution
    qcov = eye(2);
    
    %% Define passed structs
    
    clear model data params options
    
    model.ssfun = forward;
    
    % data, passed straight to forward
    CM = [5,2.5;2.5,5];
    ICM = inv(CM);
    mu = [2.5,7.5];
    data = struct('mu',mu,'ICM',ICM);
    
    params.par0 = start;
    
    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    %% Call to DRAM
    
    [results,chain] = dramrun(model,data,params,options);
    
    %% Plotting
    
    plotbivariate(chain(:,1),chain(:,2),'X','Y')
    
end