dataLength = length(fts_window.dates);
windowSize = 252;

zetas_delay5 = NaN(dataLength - windowSize,1);
zetas_mean10 = NaN(dataLength - windowSize,1);
zetas_pad20 = NaN(dataLength - windowSize,1);
zetas_win10 = NaN(dataLength - windowSize,1);
spreads = NaN(dataLength - windowSize,1);
rvols10 = NaN(dataLength - windowSize,1);
range = windowSize+1:dataLength;
for t = range
    coefs = calcwt(fts_window((t - windowSize):(t-1)));
    tresh = calthreshold(coefs,'global_mean');
    zs = ones(10,1);
    for j = 1:10
        cf = coefs(:,end-j+1);
        zs(j) = calzeta2(cf, tresh);
    end
    zetas_delay5(t - windowSize) = zs(4);
    zetas_mean10(t - windowSize) = mean(zs);
    
    lookbackprcs  = fts2mat(fts_window.series1(t-20:t)); % 21 days prices up to today
    lookbackrets = diff(log(lookbackprcs)); % 20 returns up to today
    rvol = sqrt(252/20*sum(lookbackrets.^2));
    rvols10(t - windowSize) = sqrt(252/10*sum(lookbackrets(end-9:end).^2));
    ivol = ivol_window.series1(datestr(fts_window.dates(t-20))); % fts2mat fail when query is empty.
    if ~isempty(ivol)
        spreads(t - windowSize) = 100*rvol - fts2mat(ivol);
    end
    
    % zeta with padding
    coefs = calcwt2(fts_window((t - windowSize):(t-1)), 20);
    tresh = calthreshold(coefs,'global_mean');
    zetas_pad20(t - windowSize) = calzeta2(coefs(:,end), tresh);
    
    % zeta with windowing average
    winfunc = @(x) exp(-0.5*(x-1));
    zetas_win10(t - windowSize) = zs'*winfunc(1:10)'/sum(winfunc(1:10));
end

testPeriod = fts_window.dates(range);

