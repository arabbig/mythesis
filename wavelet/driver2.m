set(0,'DefaultFigureWindowStyle','docked')

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
dataset = 'dax';
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
fprintf('>>>Windowed-zeta\n');
func_volplot2(tsobj.dates,fts2mat(tsobj.zetasPD20),fts2mat(tsobj.rv_imv));
func_spreadtab(fts2mat(tsobj.zetasPD20),fts2mat(-tsobj.rv_imv));