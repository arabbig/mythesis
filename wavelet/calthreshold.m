function thr = calthreshold(coefs,str)
meants = mean(coefs);
medts = median(coefs);
switch str
    case 'global_mean'
        thr = mean(meants);
    case 'global_sdmed'
        thr = std(medts);
    case 'local_mean'
        thr = mean(meants(end-20:end));%+std(meants(end-20:end));
    case 'local_med'
        thr = std(medts(end-20:end));
end