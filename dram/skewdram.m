function [results,chain,s2chain] = skewdram(model,data,params,options)
% DRAMRUN  Metropolis-Hastings MCMC run with adaptive delayed rejection (DRAM)
%
% This function generates MCMC chain using DRAM adaptation for a model defined
% by user supplied sum-of-squares function and with additive i.i.d. Gaussian
% errors for the observations. The error variance sigma2 is updated
% using conjugate inverse gamma distribution.
%
% This function will be modified for ABC-MCMC and ABC-DRAM.
%
% [results,chain,s2chain] = dramrun(model,data,params,options)
%
% input:
%
% model.ssfun =    ; % sum-of-squares function handle,
%                    % that returns  -2*log(p(y|par))
% model.priorfun = ; % prior "sum-of-squares" function handle
%                    % that returns -2*log(p(par)),
%                    % default: inline('0','x','params')
%
% data   = ;         % extra argument for ssfun (to pass the data etc.)
%
% params.par0   =  ; % initial parameter vector (a row vector)
% params.sigma2 =  1;% initial/prior value for the Gaussian error variance
% params.n0     = -1;% precision of sigma2 as imaginative observations
%                    %   if n0<0, no sigma2 update
% params.n      = ;  % number of actual observations (for sigma2 update)
% params.bounds = ;  % 2*npar matrix of parameter bounds
%                    % default: [-Inf,Inf]
%
% options.nSimu  = 2000;   % length of the chain
% options.qcov   = ;       % proposal covariance matrix
%
% parameters for DRAM
% options.adaptint = 10;  % how often to adapt, if zero, no adaptation
% options.drscale  = 3;   % scale for the second proposal, if zero, no DR
%
% output:
%
% results  structure that contains some info about the run
% chain    nsimu*npar MCMC chain
% s2chain  sigma� chain (if generated)

% calls covupd.m for covariance update and (optionally) gammar_mt.m for
% gamma variates

% this is a 'simple' version for demonstration and educational purposes

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.0 $  $Date: $

% Modified
% Tom Connell <thomas.connell@hdr.mq.edu.au>
% 2018


%% Required Variables (IMPORTANT)

% least-squares function - ssfun(parameters,data) - calls forward problem,
% retuns error which is fed directly into the evaluations of the acceptance
% probability
ssfun  = getpar(model,'ssfun');

% initial parameter vector (passed as array)
par0 = getpar(params,'par0'); 
par0 = par0(:)'; % row vector

% proposal covariance
qcov = getpar(options,'qcov');


%% Define MCMC variables from structs

% calls getpar(struct, field, default value)

% number of Markov chain Monte Carlo timesteps/simulations
nSimu  = getpar(options,'nsimu',10000);

% number of parameters - defined by size of par0
nPar = length(par0);

% 2*npar matrix of parameter bounds 
% (theta1 min,theta2 min; theta1 max, theta2 max) for 2*nPar
bounds = getpar(params,'bounds',(ones(nPar,2)*diag([-Inf,Inf]))');

% prior "sum-of-squares", -2*log(p(theta)) - calls prior and is fed
% directly into the acceptance probability
priorfun = getpar(model,'priorfun',@(x) 0);


%% DRAM parameters

% how often to adapt, if zero, no adaptation
adaptInt = getpar(options,'adaptint',100);

% scale for the second proposal, if zero, no DR
drScale  = getpar(options,'drscale',3);

% scale for adapting the propsal
adaScale = getpar(options,'adascale',2.4/sqrt(nPar));

% blow factor for covariace update
qcovadj  = getpar(options,'qcovadj',1e-5);


% precision of sigma2 as imaginative observations
%  if n0<0, no sigma2 update
n0  = getpar(params,'n0',-1);

% initial/prior value for the Gaussian error variance
sigma2 = getpar(params,'sigma2',1);

% number of observations (needed for sigma2 update)
if n0>=0
    n = getpar(params,'n');
end

% defines the binary values dodr, if dodr = 1 then DR is used. If dodr = 0
% then DR is not used. 
if drScale<=0
    dodr=0;
else
    dodr=1;
end

% how often to output statistics about the Markov chain
printint  = getpar(options,'printint',500);

% Cholesky factor of proposal covariance
R = chol(qcov); 

% If running DR, compute the proposal matrix for the second try
if dodr
    % second proposal for DR try
    R2 = R./drScale;
    iR = inv(R);
end

% we store the chain here, note one row per time step, NOT column
chain = zeros(nSimu,nPar);  

s20 = 0;
% If an error variance chain is being run...
if n0>=0
    % the sigma2 chain
    s2chain = zeros(nSimu,1);
    s20 = sigma2;
else
    s2chain = [];
end


%% Initialise values required for begginging the chain

% first row of the chain
oldpar = par0(:)';

% first sum-of-squares
oldss = ssfun(oldpar,data,0);

% first prior evaluation
oldprior = priorfun(oldpar);

acce = 1; % how many accepted moves

chain(1,:) = oldpar;

% If we are running an error variance chain...
if s20>0
  s2chain(1,:) = sigma2;
end

% covariance update uses these to store previous values
chaincov = []; chainmean = []; wsum = []; lasti = 0;


%% the simulation loop

for isimu=2:nSimu
    
    % information about the chain state output on every printint
    if isimu/printint == fix(isimu/printint)
        fprintf('isimu=%d, %d%% done, accepted: %d%%\n',...
            isimu,fix(isimu/nSimu*100),fix((acce/isimu)*100));
    end
    
    %% Propose new move
           
    % gets set to an index from 1 to 20
    highest_misfit_data = stochasticworstdata(oldpar,data);

    change = worstdatadirection(highest_misfit_data,oldpar,data);
    
    which_par = jcwhichpartochange(highest_misfit_data);
    
    update = zeros(1,nPar);
    clear q
    q = 0;    
    for ii = 1:nPar
        if which_par(1,ii)
            update(1,ii) = rskewnorm(0,qcov(ii,ii),change);
            % top term
            q = skewnormal(oldpar(1,ii),oldpar(1,ii)+update(1,ii),250,change);
            q = q/skewnormal(oldpar(1,ii)+update(1,ii),oldpar(1,ii),250,change);
        end
    end

    newpar = oldpar;

    %newpar(which_par) = newpar(which_par)+update(which_par);
    newpar = newpar + (update.*which_par);
    
    accept = 0;
    
    %% first M-H evaluation
    % check bounds
    if any(newpar<bounds(1,:)) || any(newpar>bounds(2,:))
        newss = Inf;
        newprior = 0;
        alpha12 = 0;
    else % inside bounds, check if accepted
        newss  = ssfun(newpar,data,0); % sum-of-squares
        newprior = priorfun(newpar); % prior ss
    
        %alpha12 = min(1,exp(-0.5*(newss-oldss)/sigma2 -0.5*(newprior-oldprior)));
        alpha12 = acceptanceRule(1,exp(-0.5*(newss-oldss)/sigma2 -0.5*(newprior-oldprior)));
        
        alpha12 = alpha12*q;
    
        if rand < alpha12 % we accept
            accept   = 1;
            acce     = acce+1;
            oldpar   = newpar;
            oldss    = newss;
            oldprior = newprior;
        end
    end
    
    % write misfit
    ssfun(oldpar,data,1); 
    
    %% DR
    if accept == 0 && dodr % we reject, but make a new try (DR)
        
        % a new try
        newpar2 = oldpar+randn(1,nPar)*R2;  
        
        % second M-H evaluation
        % check bounds
        if any(newpar2<bounds(1,:)) || any(newpar2>bounds(2,:))
            newss2 = Inf;
            newprior2 = 0;
        else % inside bounds
            
            newss2    = ssfun(newpar2,data);
            newprior2 = priorfun(newpar2);
            
            % evaluate delayed rejection acceptance probability
            alpha32 = min(1,exp(-0.5*(newss-newss2)/sigma2 -0.5*(newprior-newprior2)));
            
            l2 = exp(-0.5*(newss2-oldss)/sigma2 - 0.5*(newprior2-oldprior));
            q1 = exp(-0.5*(norm((newpar2-newpar)*iR)^2-norm((oldpar-newpar)*iR)^2));
            
            alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
            
            if rand < alpha13 % we accept
                accept = 1;
                acce = acce+1;
                oldpar = newpar2;
                oldss  = newss2;
                oldprior = newprior2;
            end
        end
    end
    
    % store current chain state
    chain(isimu,:) = oldpar;
    
    % update the error variance sigma2
    if s20 > 0
        sigma2  = 1./gammar_mt(1,1,(n0+n)./2,2./(n0*s20+oldss));
        s2chain(isimu,:) = sigma2;
    end
    
    
    %% AM
    if adaptInt>0 && fix(isimu/adaptInt) == isimu/adaptInt && isimu/adaptInt > 1
        
        % update covariance and mean of the chain
        [chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
            chaincov,chainmean,wsum);
        
        lasti = isimu;
        
        [Ra,is] = chol(chaincov + eye(nPar)*qcovadj);
        
        % singular cmat
        if is 
            fprintf('Warning cmat singular, not adapting\n');
        else
            R = Ra*adaScale;
            if dodr
                R2 = R./drScale; % second proposal for DR try
                iR = inv(R);
            end
        end
    end
    
    
end


%% Calculate output information struct

[chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
                                   chaincov,chainmean,wsum);

results.class = 'MCMC results';
results.accepted=acce./nSimu; % acceptance ratio
results.mean = chainmean;
results.cov  = chaincov;
results.qcov = R'*R;
results.R = R;
results.nsimu = nSimu;
results.drscale = drScale;
results.adascale = adaScale;
results.adaptint = adaptInt;


%% parameter setting function
function y=getpar(options,par,default)
%GETPAR get parameter value from a struct
% options   options struct
% par       parameter value to extract from the struct
% default   default value if par is not a member of the options struct

if isfield(options,par)
  y = getfield(options,par);
elseif nargin>2
  y = default;
else
  error(sprintf('Need value for option: %s',par));
end