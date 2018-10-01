function [passed_chain] = bananatestfunction(method_to_use)

% Test DRAM/MCMC with a banana shaped distribution

% run all the methods
methods = {'MH','AM', 'DR', 'DRAM'};
for mi=method_to_use:method_to_use
  
  method = methods{mi};

  nsimu = 1000002;
                                                          
  switch method
   case 'MH'
    drscale  = 0;
    adaptint = 0;
   case 'DR'
    drscale  = 2;
    adaptint = 0;
   case 'AM'
    drscale  = 0;
    adaptint = nsimu/6;
   case 'DRAM'
    drscale  = 2;
    adaptint = nsimu/6;
  end

  mu   = [0 0];         % center
  cmat = [1 0.9;0.9 1]; % target covariance
  imat = inv(cmat);

  %bpar = [1.3 1.5]; % "bananity" of the target, see bananafun.m
  bpar = [1 1]; % "bananity" of the target, see bananafun.m
  % "sum-of-squares" function as inline function
  bananass = @(x,d) bananafun(x-d.mu,d.bpar,0)*d.imat*bananafun(x-d.mu,d.bpar,0)';
  %bananass = inline('bananafun(x-d.mu,d.bpar,1)*d.imat*bananafun(x-d.mu,d.bpar,1)''','x','d');

  % create input arguments for the dramrun function
  clear model params options

  model.ssfun    = bananass;

  params.par0    = mu; % initial value

  data = struct('mu',mu,'imat',imat,'bpar',bpar);   

  options.nsimu    = nsimu;
  options.adaptint = adaptint;
  options.drscale  = drscale;
  options.qcov     = eye(2)*5; % initial covariance 

  [~,chain] = dramrun(model,data,params,options);
  
end

passed_chain = chain;