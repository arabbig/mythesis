%% Calculating CWT, Zetas
N = min(length(fts_window));
windowSize = 252;

zetas_delay5 = ones(N-windowSize, 1);
zetas_mean10 = ones(N-windowSize, 1);
%zetas_max5 = ones(N-windowSize, 1);
%zetas_max10 = ones(N-windowSize, 1);
RANGE = (windowSize + 1) : N ;
for date = RANGE 
    coefs = calcwt(fts_window((date - windowSize):(date-1)));
    tresh = calthreshold(coefs,'global_mean');
    zs = ones(10,1);
    for j = 1:10
        cf = coefs(:,end-j+1);
        zs(j) = calzeta2(cf, tresh); 
    end
    zetas_delay5(date-windowSize) = zs(4);
    zetas_mean10(date-windowSize) = mean(zs);    
end