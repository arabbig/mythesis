%%
load zetaLeadSpread.mat
dataset = 'spx';
tsobj = outputTS.(dataset);
set(0,'DefaultFigureWindowStyle','docked')
%%
% Delay5-zeta 
% [+] It has longer days predictive horizons than ten days vol (esp. at lag 15)
fprintf('>>>Delay5-zeta\n');
func_volplot2(tsobj.dates,fts2mat(tsobj.zetasDL5),fts2mat(tsobj.rv_imv));
figure;
func_spreadtab(fts2mat(tsobj.zetasDL5),fts2mat(-tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasDL5),fts2mat(-tsobj.rv_imv),15);
crosscorr(fts2mat(tsobj.rv10),fts2mat(tsobj.rv_imv));
title('xcf(RV_{10},RV_{20})','FontSize',18);
% Mean10-zeta 
% [-] It has shorter predictive horizons
% [+] Very nice time series
fprintf('>>>Mean10-zeta\n');
func_volplot2(tsobj.dates,fts2mat(tsobj.zetasMN10),fts2mat(tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasMN10),fts2mat(-tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasMN10),fts2mat(-tsobj.rv_imv),15);

% Windowed-zeta
% fprintf('>>>Windowed-zeta\n');
% func_volplot2(tsobj.dates,fts2mat(tsobj.zetasWN10),fts2mat(tsobj.rv_imv));
% func_spreadtab(fts2mat(tsobj.zetasWN10),fts2mat(-tsobj.rv_imv));

% Lastday-zeta
% fprintf('>>>Nodelay_zeta\n');
% func_volplot2(tsobj.dates,fts2mat(tsobj.zetasPD20),fts2mat(tsobj.rv_imv));
% func_spreadtab(fts2mat(tsobj.zetasPD20),fts2mat(-tsobj.rv_imv));

%% show in writeup_midjune
t = 150;
span = 500;

figure;
subplot(2,1,1);
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t:t+span)),tsobj.dates(t:t+span),fts2mat(tsobj.index(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t) and Price','FontSize',14);

subplot(2,1,2);
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t+5:t+5+span)),tsobj.dates(t:t+span),fts2mat(tsobj.index(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t+5) and Price','FontSize',14);



