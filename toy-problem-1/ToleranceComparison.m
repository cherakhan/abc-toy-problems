% A comparison of the impact changing the tolerance has on an ABC rejection
% algorithm. This is a repition of the analysis in Gaussian1dRejection.m
% for varying tolerance. The algorithm targets mu and sigma in a Gaussian
% distibution. The summary stats are sample mean and sample s.t.d

clear;clc;clf


%% Generate "observed data"

observed_data = (randn(100,1)*2)+5;
observed_sumstats = [mean(observed_data); std(observed_data)];


%% Simulate from prior
% p(mu) = U(0,10)
% p(sigma^2) = U(0,100)

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

for jj = 1:3
    
    clear posteriorMu posteriorSigma
        
    if jj == 1
        tolerance = 0.1;
    end
    if jj == 2
        tolerance = 0.5;
    end
    if jj == 3
        tolerance = 1;
    end


    %% Perform rejection

    postmu = zeros(n_sims,1);
    postsigma = zeros(n_sims,1);

    for ii = 1:n_sims
       if abs(sim_sumstats(ii,1) - observed_sumstats(1))...
               + abs(sim_sumstats(ii,2) - observed_sumstats(2)) <= tolerance

           postmu(ii) = prior_mu(ii);
           postsigma(ii) = prior_sigma(ii);
       end
    end

    posterior_mu = postmu(postmu~=0);
    posterior_sigma = postsigma(postsigma~=0);
    
    if jj == 1
        [bandwidth1,MuDensityP1,MuMeshP1] = kde(posterior_mu,2^7,[0,10]);
        [bandwidth2,SigmaDensityP1,SigmaMeshP1] = kde(posterior_sigma,2^7,[0,10]);
    end
    if jj == 2
        [bandwidth1,MuDensityP5,MuMeshP5] = kde(posterior_mu,2^7,[0,10]);
        [bandwidth2,SigmaDensityP5,SigmaMeshP5] = kde(posterior_sigma,2^7,[0,10]);
    end
    if jj == 3
        [bandwidth1,MuDensity1,MuMesh1] = kde(posterior_mu,2^7,[0,10]);
        [bandwidth2,SigmaDensity1,SigmaMesh1] = kde(posterior_sigma,2^7,[0,10]);
    end
    
end

%% Plotting
figure(1)

colors = brewermap(3,'Set1');

subplot(1,2,1)
plot(MuMeshP1,MuDensityP1,'color',colors(1,:),'linewidth',1)
hold on
plot(MuMeshP5,MuDensityP5,'color',colors(2,:),'linewidth',1)
plot(MuMesh1,MuDensity1,'color',colors(3,:),'linewidth',1)
xlabel('\mu')
ylabel('posterior density')
xlim([0,10])
text(0.05,0.95,'a','Units','normalized')
set(gca,'Ytick',[])
xlim([3.25,6.75])

subplot(1,2,2)
plot(SigmaMeshP1,SigmaDensityP1,'color',colors(1,:),'linewidth',1)
hold on
plot(SigmaMeshP5,SigmaDensityP5,'color',colors(2,:),'linewidth',1)
plot(SigmaMesh1,SigmaDensity1,'color',colors(3,:),'linewidth',1)
xlabel('\sigma')
ylabel('posterior density')
xlim([0,10])
legend({'\epsilon = 0.1','\epsilon = 0.5','\epsilon = 1'})
legend boxoff
text(0.05,0.95,'b','Units','normalized')
set(gca,'Ytick',[])
xlim([0.5,3.25])

set(gcf,'units','centimeters','position',[0,0,25,10],'papersize',[25,10])
print -dpdf -painters tolerancecomparison.pdf