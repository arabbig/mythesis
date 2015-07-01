%% Plotting SW-AR results
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/chanp/Documents/MATLAB/brml/')
setup
load ../sw_ar_2001_2015.mat
load('../wavelet_2001_2015.mat')
dataset = 'spx';
tsobj = combineTS.(dataset);
close = fts2mat(tsobj.close);
close(1) =[];
clear combineTS;
%%
figure;
subplot(2,1,1);
[axs,h1,h2] = plotyy(dates,abs(ret),dates,rv_imv);
datetick(axs(1));
datetick(axs(2));
line(repmat(dates(1755),1,2),[-100, 100],'Color','black','Parent',axs(2))
set(axs(1),'XLim',[min(dates),max(dates)]);
set(axs(2),'XLim',[min(dates),max(dates)]);
set(axs(2),'YLim',[-100,50]);
set(axs(1),'YLim',[-0.10,0.35]);
%set([h1;h2],'Marker','.')
set(axs(1),'ycolor','r')
set(axs(1),'ycolor','b')  % Left color red, right color blue...

colormap bone;
subplot(2,1,2);
tmp = [ret_prob, brml.condexp(ret_logalpha)];
imagesc(tmp);
[~,ret_smax] = max(tmp); 

m1 = mean(abs(ret(ret_smax==1)));
m2 = mean(abs(ret(ret_smax==2)));
m3 = mean(abs(ret(ret_smax==3)));
[~,srt] = sort([m1,m2,m3]);

col = zeros(1, length(ret));
mycmap = [1,5,10];
for ll = 1:3
    col(ret_smax==srt(ll)) = mycmap(ll);
end

figure;
subplot(3,1,1);
scatter(dates,ret,[],col);
datetick();
subplot(3,1,2);
scatter(dates,close,[],col);
datetick();
subplot(3,1,3);
scatter(dates,rv_imv,[],col);
datetick();


% figure;
% subplot(2,1,1);
% scatter(dates(ret_smax==srt(2)),rv_imv(ret_smax==srt(2)),'MarkerEdgeColor',[0.3 0.5 1]);
% hold on
% scatter(dates(ret_smax==srt(1)),rv_imv(ret_smax==srt(1)),'MarkerEdgeColor',[0.4 0.9 0.6]);
% hold on
% scatter(dates(ret_smax==srt(3)),rv_imv(ret_smax==srt(3)),'MarkerEdgeColor',[1 0.5 0.5]);
% datetick();
% subplot(2,1,2);
% scatter(dates(ret_smax==srt(2)),(ret(ret_smax==srt(2))),'MarkerEdgeColor',[0.3 0.5 1]);
% hold on
% scatter(dates(ret_smax==srt(1)),(ret(ret_smax==srt(1))),'MarkerEdgeColor',[0.4 0.9 0.6]);
% hold on
% scatter(dates(ret_smax==srt(3)),(ret(ret_smax==srt(3))),'MarkerEdgeColor',[1 0.5 0.5]);
% datetick();

[~,~,rnk]=unique([m1,m2,m3]);
str = {'low','med','hig'};
disp('P(s(t) = i| s(t-1) = j)')
fprintf('    %s      %s      %s\n',char(str(rnk(1))),char(str(rnk(2))),char(str(rnk(3))));
fprintf('%s %f %f %f\n',char(str(rnk(1))),ret_stran(1,1),ret_stran(1,2),ret_stran(1,3));
fprintf('%s %f %f %f\n',char(str(rnk(2))),ret_stran(2,1),ret_stran(2,2),ret_stran(2,3));
fprintf('%s %f %f %f\n',char(str(rnk(3))),ret_stran(3,1),ret_stran(3,2),ret_stran(3,3));

figure;

subplot(2,1,1);
scatter(dates(ret_smax==srt(2)),(close(ret_smax==srt(2))),'MarkerEdgeColor',[0.3 0.5 1]);
hold on
scatter(dates(ret_smax==srt(1)),(close(ret_smax==srt(1))),'MarkerEdgeColor',[0.4 0.9 0.6]);
hold on
scatter(dates(ret_smax==srt(3)),(close(ret_smax==srt(3))),'MarkerEdgeColor',[1 0.5 0.5]);
hold on
datetick();


subplot(2,1,2);
scatter(dates(ret_smax==srt(2)),(ret(ret_smax==srt(2))),'MarkerEdgeColor',[0.3 0.5 1]);
hold on
scatter(dates(ret_smax==srt(1)),(ret(ret_smax==srt(1))),'MarkerEdgeColor',[0.4 0.9 0.6]);
hold on
scatter(dates(ret_smax==srt(3)),(ret(ret_smax==srt(3))),'MarkerEdgeColor',[1 0.5 0.5]);
hold on
datetick();


%%
figure;
index = 1:length(ret);
[~,rlow_smax] = max(rlow_prob); 

x1 = mean(abs(ret(rlow_smax==1)));
x2 = mean(abs(ret(rlow_smax==2)));
[~,lev] = sort([x1,x2]);

rlow_probfwd = brml.condexp(rlow_logalpha);
col = zeros(1,length(ret));
col(1756:end) = rlow_probfwd(lev(2),:);
subplot(4,1,1);
scatter(dates,close,[],col);
datetick();
subplot(4,1,2);
scatter(dates,ret,[],col);
datetick();
subplot(4,1,3);
scatter(dates,rv_imv,[],col);
datetick();
subplot(4,1,4);
scatter(dates,col,[],col);
datetick();

%% Train on spread itself
% unsuccesful attemp
% shown here for completeness
% 
% figure;
% subplot(2,1,1);
% [axs,h1,h2] = plotyy(dates,abs(ret),dates,rv_imv);
% datetick(axs(1));
% datetick(axs(2));
% line(repmat(dates(1755),1,2),[-100, 100],'Color','black','Parent',axs(2))
% set(axs(1),'XLim',[min(dates),max(dates)]);
% set(axs(2),'XLim',[min(dates),max(dates)]);
% set(axs(2),'YLim',[-100,50]);
% set(axs(1),'YLim',[-0.10,0.35]);
% %set([h1;h2],'Marker','.')
% set(axs(1),'ycolor','r')
% set(axs(1),'ycolor','b')  % Left color red, right color blue...
% 
% colormap bone;
% subplot(2,1,2);
% tmp = [spread_prob, brml.condexp(spread_logalpha)];
% imagesc(tmp);
% [~,spread_smax] = max(tmp); 
% 
% m1 = mean((rv_imv(spread_smax==1)));
% m2 = mean((rv_imv(spread_smax==2)));
% m3 = mean((rv_imv(spread_smax==3)));
% [~,srt] = sort([m1,m2,m3]);
% figure;
% scatter(dates(spread_smax==srt(2)),rv_imv(spread_smax==srt(2)),'MarkerEdgeColor',[0.3 0.5 1]);
% hold on
% scatter(dates(spread_smax==srt(1)),rv_imv(spread_smax==srt(1)),'MarkerEdgeColor',[0.4 0.9 0.6]);
% hold on
% scatter(dates(spread_smax==srt(3)),rv_imv(spread_smax==srt(3)),'MarkerEdgeColor',[1 0.5 0.5]);
% datetick();
