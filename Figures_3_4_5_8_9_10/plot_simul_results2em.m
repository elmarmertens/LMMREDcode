burn_in=1000;

yirf_simul_full=yirf_simul;
yirf_simul=yirf_simul(:,burn_in+1:end,:);

[Nyy, ~, Nsol] = size(yirf_simul);
Nlag = 20;

mean_simul=squeeze(mean(yirf_simul,2));
std_simul=squeeze(std(yirf_simul,0,2));
std_simul_full_info=squeeze(std(yirf_simul_full_info,0,2));

% allocate memory
simul_ac           = NaN(Nyy, Nsol, Nlag + 1);
simul_ac_full_info = NaN(Nyy, Nlag + 1);

for kk=1:Nyy
    for jj=1:Nsol
        simul_ac(kk,jj,:)=autocorr(yirf_simul(kk,:,jj), Nlag);
    end
end

for kk=1:Nyy
    simul_ac_full_info(kk,:)=autocorr(yirf_simul_full_info(kk,:), Nlag);
end



%% plot: top row 
figure;
for jj=1:3
    subplot(2,3,jj)
    for ee=1:size(yirf_simul,3)
        plot(squeeze(simul_ac(jj,ee,2:end)),'r','LineWidth',1.5)
        grid on
        hold on
        axis tight
    end
    plot(squeeze(simul_ac_full_info(jj,2:end)),'b','LineWidth',2)
    xlabel('horizon')
    title([' autocorrelation of ', YYlabels{jj}])
end

%% bottom panels

minni = NaN(3,1);
maxi  = NaN(3,1);

subplot(2,3,4:6)
hold on
for jj=1:3
    
    minni(jj) = min(std_simul(jj,:)/std_simul_full_info(jj));
    maxi(jj)  = max(std_simul(jj,:)/std_simul_full_info(jj));
    
    plot([jj jj], [minni(jj) maxi(jj)], 'r-', 'LineWidth',20)
    
    
end
set(gca, 'xtick', 1 : 3)
set(gca, 'xticklabel', YYlabels(1:3))
xlim([0 4])

plot(xlim, [1 1], 'k:')
title('rel std')
% legend(YYlabels{1:3})

figname = 'figure_10';
print(figname, '-depsc')
savefig(figname)

