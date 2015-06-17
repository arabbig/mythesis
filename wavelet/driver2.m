set(0,'DefaultFigureWindowStyle','docked')
%%
for dataindex = 1:4
    main_readdata;
    main_calzeta2; 
    range = windowSize+1:dataLength;
    indexTS = chfield(fts_window(range), 'series1', 'index');
    zetaTS = fints(testPeriod,[zetas_delay5, zetas_mean10, zetas_pad20, zetas_win10], {'zetasDL5','zetasMN10','zetasPD20','zetasWN10'});
    spreadTS = fints(testPeriod, spreads, 'rv_imv');
    rvols10TS = fints(testPeriod, rvols10, 'rv10');
    outputTS.(dataset) = merge(indexTS,zetaTS,spreadTS,rvols10TS);
end

save zetaLeadSpread.mat outputTS

%%
dataset = 'spx';
tsobj = outputTS.(dataset);

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
fprintf('>>>Windowed-zeta\n');
func_volplot2(tsobj.dates,fts2mat(tsobj.zetasWN10),fts2mat(tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasWN10),fts2mat(-tsobj.rv_imv));

% Lastday-zeta
fprintf('>>>Nodelay_zeta\n');
func_volplot2(tsobj.dates,fts2mat(tsobj.zetasPD20),fts2mat(tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasPD20),fts2mat(-tsobj.rv_imv));

%% show in writeup_midjune
t = 150;
span = 200;

figure;
subplot(3,1,1);
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t:t+span)),tsobj.dates(t:t+span),fts2mat(tsobj.index(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t) and Price','FontSize',14);

subplot(3,1,2);
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t+5:t+5+span)),tsobj.dates(t:t+span),fts2mat(tsobj.index(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t+5) and Price','FontSize',14);

subplot(3,1,3);
t = 200;
span = 500;
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t:t+span)),tsobj.dates(t:t+span),fts2mat(tsobj.rv_imv(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
set(axs(2),'YLim',[-20,20]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t) and vol-spread','FontSize',14);
hl = legend([h1;h2],'\zeta','\sigma_{real}-\sigma_{implied}');

%% segmentation of price
t = 1;
span = 500;
[pks,locs]=findpeaks(fts2mat(tsobj.zetasDL5(t:t+span))); % locs ranges in 1:1+span.

subplot(2,1,1)
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t+5:t+5+span)),tsobj.dates(t:t+span),fts2mat(tsobj.index(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t+5) and Price','FontSize',14);

for ii = 1:length(locs)
    pos = locs(ii) + t - 1; % convert to index of time series
    if pos>=11 && pos+10 <= t+span && pks(ii) >= max(fts2mat(tsobj.zetasDL5(pos-10:pos+10)))
        line([tsobj.dates(pos-5) tsobj.dates(pos-5)],ylim,'Color','green','LineStyle',':');
        line([tsobj.dates(pos-6) tsobj.dates(pos-6)],ylim,'Color','magenta','LineStyle',':');
    end
end

subplot(2,1,2);
[axs,h1,h2]=plotyy(tsobj.dates(t:t+span),fts2mat(tsobj.zetasDL5(t+5:t+5+span)),tsobj.dates(t:t+span),fts2mat(tsobj.rv_imv(t:t+span)));
hold on;
set(axs(1),'YLim',[0,1.5]);
set(axs(2),'YLim',[-20,20]);
datetick(axs(1));
datetick(axs(2));
title('\zeta(t) and vol-spread','FontSize',14);
hl = legend([h1;h2],'\zeta','\sigma_{real}-\sigma_{implied}');
refline(axs(2),[0 3]);
% for ii = 1:length(locs)
%     pos = locs(ii) + t - 1; % convert to index of time series
%     if pos>=11 && pos+10 <= t+span && pks(ii) >= max(fts2mat(tsobj.zetasDL5(pos-10:pos+10)))
%         line([tsobj.dates(pos) tsobj.dates(pos)],ylim,'Color','green');
%         line([tsobj.dates(pos-5) tsobj.dates(pos-5)],ylim,'Color','magenta');
%     end
% end

