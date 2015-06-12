function func_spreadtab(zetas,iv_rv,lag)
    X = [zetas -iv_rv];
    X(any(isnan(X),2),:) = []; % remove missing data
    nonmissing_zetas  = X(:,1);
    nonmissing_rv_iv  = X(:,2);

    if nargin < 3
        % find best lag. use rv_iv instead of iv_rv
        [xcf,lags] = crosscorr(nonmissing_zetas, nonmissing_rv_iv);
        xcf(lags < 0) = -Inf;
        [~,mi] =max(xcf); 
        lag = lags(mi);
    end
    
    
    zetasLag = zetas(1:end-lag);
    iv_rv(1:lag) = [];

    fprintf('E(impvol-rvol | z(t-l) > X), l = %d\n', lag);
    fprintf('Days during testing = %d\n',length(zetasLag));
    fprintf('   X      days     AvgPayoff\n');
    for X = linspace(0,0.9,10)
        ii  = zetasLag > X;
        fprintf('%6.3f%9.d%11.3f\n', X, sum(ii)+0,nanmean(iv_rv(ii)));
    end
end
