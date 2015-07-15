set(0,'DefaultFigureWindowStyle','docked')

%% Segmentation of price

% This chunk is prototype,i.e., it does operate on 2008+ data

load('../wavelet_2001_2015.mat')
load ../pricedata_fullperiod.mat;

addpath('/Users/chanp/Documents/MATLAB/brml/')
setup

dataset = {'spx','ukx','dax','nky'};

for di = 2:2
    
tsobj = combineTS.(char(dataset(di)));

%% Dividing data

pfts = asset.(char(dataset(di))).priceData(datestr(tsobj.dates));
ifts = asset.(char(dataset(di))).ivolData(datestr(tsobj.dates));

ret = pfts.ret;
iv = chfield(ifts.ivol_percent,'ivol_percent','iv');
rv = chfield(100*pfts.rvar_lead.^0.5','rvar_lead','rv');
both = merge(iv, rv,'DateSetMethod','Intersection');
both.payoff = fts2mat(both.iv) - fts2mat(both.rv);
payoff = both.payoff;

%ret = chfield(diff(log(tsobj.close)),'close','ret');
%payoff = -chfield(leadts(tsobj.rv_imv,20),'rv_imv','payoff');

tsobj = merge(tsobj,ret,payoff);
tsobj = tsobj(2:end-20); % the first day is NaN for ret, last 20 days you don't have the payoff

trainTS = tsobj('::11/30/2007');
train_ret = fts2mat(trainTS.ret);

testTS = tsobj('12/01/2007::');
test_ret = fts2mat(testTS.ret);
test_payoff = fts2mat(testTS.payoff);

%% Train 3-state HMM

L = 1; S = 3; Tskip = 5;
opts.plotprogress = false;
opts.maxit = 20;

% datestr(dates(1755)) is 31 Dec 2008
[ret_ar,ret_sigma2,ret_stran,ret_prob]=brml.SARlearn((train_ret)', L, S, Tskip, opts);
ret_prior = brml.condp(sum(ret_prob,2));

% Estimation of the most likely state on each date
[~, smax1] = max(ret_prob); 

m1 = mean((train_ret(smax1==1)).^2);
m2 = mean((train_ret(smax1==2)).^2);
m3 = mean((train_ret(smax1==3)).^2);
[~,lev1] = sort([m1,m2,m3]);


%% Train 2-state HMM

% rlow = train_ret(smax1==lev1(1));
% L = 1; S = 2; Tskip = 5;
% opts.plotprogress = false;
% opts.maxit = 20;
% 
% [rlow_ar,rlow_sigma2,rlow_stran,rlow_prob]=brml.SARlearn(rlow', L, S, Tskip, opts);
% rlow_prior = brml.condp(sum(rlow_prob,2));

%% Making prediction (using 2-state HMM)
% [rlow_logalpha,~]=brml.HMMforwardSAR((test_ret)', rlow_stran, rlow_prior, rlow_ar, rlow_sigma2, Tskip);
% rlow_probfwd = brml.condexp(rlow_logalpha);
% 
% [~, smax2] = max(rlow_prob ); 
% x1 = mean((rlow(smax2==1)).^2);
% x2 = mean((rlow(smax2==2)).^2);
% [~,lev2] = sort([x1,x2]);
% 
% prob_highvol = rlow_probfwd(lev2(2),:);

%% Making prediction (using 3-state HMM)

[ret_logalpha,~]=brml.HMMforwardSAR((test_ret)', ret_stran, ret_prior, ret_ar, ret_sigma2, Tskip);
ret_probfwd = brml.condexp(ret_logalpha);

[~,smax3] = max(ret_probfwd);
%prob_highvol = 1-ret_probfwd(lev1(1),:);
prob_highvol = ret_probfwd(lev1(3),:);
% prob_highvol = (smax3==lev1(1)).*(1-ret_probfwd(lev1(1),:))+(smax3 ~= lev1(1)).*ret_probfwd(lev1(2),:);

%%

L = length(prob_highvol);
tailProb = zeros(L, 1);
runmax = ones(L, 1);
jump = zeros(L,1);
CALM = 0;
EXCITED = 1;
WILD = 2;
PEAK = 3;
FADE = 4;
flag = CALM;
cntWARN = 0;

for i=21:length(jump)
    v1 = prob_highvol(i-19:i);
    v1 = log(v1./(1-v1));
    mu1 = mean(v1);
    sig1 = std(v1);
    v2 = prob_highvol(i-4:i);
    tailProb(i) = normcdf((v1(end)-mu1)/sig1);
    switch flag
        case CALM
            jump(i) = 0;
            if tailProb(i) > 0.95
                flag = WILD;
            end
        case WILD
            jump(i) = 1;
            if tailProb(i) < 0.2
                flag = CALM;
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%% visual inspection
% 
% figure;
% subplot(4,1,1);
% scatter(testTS.dates,test_ret,[],jump);
% set(gca,'XMinorTick','on');
% datetick();
% xlim([testTS.dates(1), testTS.dates(end)]);
% subplot(4,1,2);
% scatter(testTS.dates, test_payoff ,[],jump);
% refline([0 0])
% set(gca,'XMinorTick','on');
% datetick();
% title('IV_{t}-RV_{t+20}');
% xlim([testTS.dates(1), testTS.dates(end)]);
% subplot(4,1,3);
% scatter(testTS.dates,prob_highvol,[],jump);
% datetick();
% xlim([testTS.dates(1), testTS.dates(end)]);
% subplot(4,1,4);
% plot(testTS.dates,tailProb);
% datetick();
% xlim([testTS.dates(1), testTS.dates(end)]);

signal.(char(dataset(di))) = fints(testTS.dates, 1 - jump, 'regimeJMP');
signal.(char(dataset(di))).zetaDL5 = tsobj.zetDL5;
signal.(char(dataset(di))).zetaPAD = tsobj.zetPAD;

end

%% save ../twosignal_2008_2015.mat signal

