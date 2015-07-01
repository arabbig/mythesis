function func_volplot2(dates,zetas,spreads)
    X = [zetas spreads];
    X(any(isnan(X),2),:) = []; % remove missing data
    nonmissing_zetas  = X(:,1);
    nonmissing_spreads  = X(:,2);

    [xcf,lags] = crosscorr(nonmissing_zetas,nonmissing_spreads);
    xcf(lags < 0) = -Inf;
    [~,mi] =max(xcf);    
    lag = lags(mi);
    
    t = 21;
    span = 1500;
    
    figure;
    subplot(2,1,1);
    crosscorr(nonmissing_zetas,nonmissing_spreads);
    title('xcf(\zeta,\sigma_{real}-\sigma_{implied})','FontSize',14,'FontWeight','bold');
    
    subplot(2,1,2);
    [hxs,h1,h2] = plotyy(dates(t:t+span), zetas(t-lag:t-lag+span),dates(t:t+span), spreads(t:t+span));
    datetick(hxs(1));
    datetick(hxs(2));
    set(hxs,'FontSize',14);
    set(hxs(1),'YLim',[0,1.5]);
    set(hxs(2),'YLim',[-20,20]);
    %set([h1;h2],'Marker','.');
    hl = legend([h1;h2],'\zeta','\sigma_{real}-\sigma_{implied}');
    set(hl,'FontSize',14);
    title(strcat('Zeta is lagged by = ',num2str(lag)));
end