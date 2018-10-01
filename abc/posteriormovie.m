function [] = posteriormovie(posterior,outputpath)

obs = csvread('bananaData.dat');

for ii = 1:length(posterior)
    figure(1)
    scatter(obs(:,1),obs(:,2),'k.')
    hold on
    bananaplot(0,0,1,1,0.9,1,1,'k-')
    bananaplot(posterior(ii,1),posterior(ii,2),posterior(ii,3),posterior(ii,4),posterior(ii,5),posterior(ii,6),posterior(ii,7),'r-')
    xlim([-10,10])
    ylim([-5,15])
    saveas(gcf,[outputpath,filesep,'post',num2str(ii),'.png'])
    clf
end