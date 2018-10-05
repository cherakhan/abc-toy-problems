% An example of an ABC rejection algorithm 
% Here we target the parameters mu and sigma for a Gaussian distribution
% The sample mean and sample standard deviation are used as summary stats

clear;clc;clf


%% Generate "observed data"

observed_data = (randn(100,1)*2)+5;
observed_sumstats = [mean(observed_data); std(observed_data)];


%% Simulate from prior
% p(mu) = U(0,10)
% p(sigma) = U(0,10)

% number of simulations
n_sims = 1000000;

prior_mu = rand(n_sims,1)*10;
prior_sigma = rand(n_sims,1)*10;


%% Forward simulate and store summaries

sim_sumstats = zeros(n_sims,2);

for ii = 1:n_sims
   
   simulated_data =  (randn(100,1)*prior_sigma(ii))+prior_mu(ii);
   
   sim_sumstats(ii,1) = mean(simulated_data);
   sim_sumstats(ii,2) = std(simulated_data); 
   
end


%% Perform rejection step

tolerance = 0.1;

% containers for posterior samples
postmu = zeros(n_sims,1);
postsigma = zeros(n_sims,1);

% metric over summary statistics
metric = @(ss1, ss2) abs(ss1 - observed_sumstats(1))...
    + abs(ss2 - observed_sumstats(2));

for ii = 1:n_sims
   
   % if the distance between summary statistics is less than the tolerance
   if metric(sim_sumstats(ii,1), sim_sumstats(ii,2)) <= tolerance
       
       % add the sample to the set of posterior samples
       postmu(ii) = prior_mu(ii);
       postsigma(ii) = prior_sigma(ii);
       
   end
   
end

% remove empty ([0,0]) entries from posterior containers
posterior_mu = postmu(postmu~=0);
posterior_sigma = postsigma(postsigma~=0);


%% Compute analytical likelihood
% Compute the likelihood/posterior analytically to compare to abc results

likelihood_calc = @(l2,var,n) ((2*pi*var)^(-n/2))*exp(-(1/(2*var))*l2);

mu = 0:0.01:10;
sigma = 0.01:0.01:10;

likelihood_table = zeros(length(mu),length(sigma));

for ii = 1:length(mu)
    for jj = 1:length(sigma)
        
        error = sum((observed_data-mu(ii)).^2);
        
        likelihood_table(ii,jj) = likelihood_calc(error,sigma(jj)^2,100);
        
    end
end

marginal_mu = zeros(length(mu),1);
for ii = 1:length(mu)
   marginal_mu(ii) = sum(likelihood_table(ii,:),2);
end

marginal_sigma = zeros(length(sigma),1);
for ii = 1:length(sigma)
   marginal_sigma(ii) = sum(likelihood_table(:,ii)); 
end

marginal_mu = marginal_mu/max(marginal_mu);
marginal_sigma = marginal_sigma/max(marginal_sigma);


%% Plot

subplot(1,3,1)
h1 = histogram(posterior_mu,'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
plot(mu,marginal_mu*max(h1.Values),'r')
set(gca,'YTick',[])
xlabel('\mu')
ylabel('posterior density')
text(0.05,0.95,'a','Units','normalized')
xlim([4,6.25])

subplot(1,3,2)
h2 = histogram(posterior_sigma,'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
plot(sigma,marginal_sigma*max(h2.Values),'r')
legend({'ABC Posterior', 'Analytical Posterior'})
legend boxoff
xlabel('\sigma')
ylabel('posterior density')
set(gca,'YTick',[])
text(0.05,0.95,'b','Units','normalized')
xlim([1.25,3])

subplot(1,3,3)
% True distribution vs. ABC posterior mean
range = 0:0.01:10;
trueDist = makedist('Normal',5,2);
abcPosterior = makedist('Normal',mean(posterior_mu),mean(posterior_sigma));

true = pdf(trueDist,range);
abc = pdf(abcPosterior,range);

plot(range,true,'b')
hold on
plot(range,abc,'k')
legend({'True model','ABC posterior mean'})
legend boxoff
set(gca,'YTick',[])
xlabel('X')
ylabel('probability density')
ylim([0,0.3])
text(0.05,0.95,'c','Units','normalized')

set(gcf,'units','centimeters','position',[0,0,30,10],'papersize',[30,10])
print -dpdf -painters rejectionsampler.pdf