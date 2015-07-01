set(0,'DefaultFigureWindowStyle','docked')

%% Segmentation of price

% This chunk is prototype,i.e., it does operate on 2008+ data

load('../wavelet_2001_2015.mat')
dataset = 'spx';
tsobj = combineTS.(dataset);
addpath('/Users/chanp/Documents/MATLAB/brml/')
setup
%% Fitting Switching AR model
ret = diff(log(fts2mat(tsobj.close)));
rv_imv = fts2mat(tsobj.rv_imv);
rv_imv(1) = [];
zet = fts2mat(tsobj.zetDL5);
zet(1) = [];
dates = tsobj.dates;
dates(1) = [];
close = fts2mat(tsobj.close);
close(1) = [];
%%

L = 1; S = 3; Tskip = 5;
opts.plotprogress = false;
opts.maxit = 20;

% v = abs(ret'); 
% a=brml.condp(randn(L,S)); % set the AR coefficients
% stran=brml.condp(ones(S,S)); % switch transition
% sprior=brml.condp(ones(S,1)); % switch prior
% sigma2=var(v)*ones(1,S);
% T=size(v,2);
% datestr(dates(1755)) is 31 Dec 2008
[ret_ar,ret_sigma2,ret_stran,ret_prob]=brml.SARlearn((ret(1:1755)'), L, S, Tskip, opts);
ret_prior = brml.condp(sum(ret_prob,2));
[ret_logalpha,~]=brml.HMMforwardSAR((ret(1756:end)'), ret_stran, ret_prior, ret_ar, ret_sigma2, Tskip);

%[spread_ar,spread_sigma2,spread_stran,spread_prob]=brml.SARlearn(rv_imv(1:1755)',L,S,Tskip,opts);
%spread_prior = brml.condp(sum(spread_prob,2));
%[spread_logalpha,~]=brml.HMMforwardSAR((rv_imv(1756:end)'), spread_stran, spread_prior, spread_ar, spread_sigma2, Tskip);


%% Refit within lowest-vol state

[~,bwd_smax] = max(ret_prob); 

m1 = mean(abs(ret(bwd_smax==1)));
m2 = mean(abs(ret(bwd_smax==2)));
m3 = mean(abs(ret(bwd_smax==3)));
[~,srt] = sort([m1,m2,m3]);

rlow = ret(bwd_smax==srt(1));

%%
L = 1; S = 2; Tskip = 5;
opts.plotprogress = false;
opts.maxit = 20;

[rlow_ar,rlow_sigma2,rlow_stran,rlow_prob]=brml.SARlearn(rlow', L, S, Tskip, opts);
rlow_prior = brml.condp(sum(rlow_prob,2));
[rlow_logalpha,~]=brml.HMMforwardSAR((ret(1756:end)'), rlow_stran, rlow_prior, rlow_ar, rlow_sigma2, Tskip);

%%
save ../sw_ar_2001_2015.mat ...
    ret rv_imv close dates rlow ...
    ret_logalpha ret_prob ret_ar ret_stran ret_sigma2 ...
    spread_logalpha spread_prob spread_ar spread_stran spread_sigma2 ...
    rlow_logalpha rlow_prob rlow_ar rlow_sigma2 rlow_stran;