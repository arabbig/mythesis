function coefs_zeta=calcwt(fts_window)


y = fts2mat(fts_window.series1);
%prev = y(1);
%y(1) = [];
%logret = [log(y(1)) - log(prev); diff(log(y))];

LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
index = y;

%%
trendType = {'poly1', 'exp1'};

f=fit(x,y, trendType{2});
figure
plot(f, x, y);
y = y - f(x);

%%
f = fit(x,y,'sin2','Upper',[Inf 1 Inf Inf 1 Inf],'Robust','on','Display','off');
figure
plot(f, x, y);
y = y - f(x);



%%
figure
subplot(3,1,1)
coefs = cwt(y,1:32,'meyr','plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');


%%

coefs_by_day = num2cell(coefs,1);
pfunc = @(q) @(x) prctile(x,q);
ninethyP = pfunc(90);
tenthP = pfunc(10);

tenthP_ccfs = cellfun(tenthP, coefs_by_day);
median_ccfs = cellfun(@median, coefs_by_day);
ninethyP_ccfs = cellfun(ninethyP, coefs_by_day);

morethan = @(tresh) @(x) min(sum(x>tresh)/length(x),0.5);
zeta = morethan(std(median_ccfs));
%morethan = @(tresh) @(x) min(sum(x>tresh)/length(x),1);
%zeta=morethan(mean(reshape(coefs,[prod(size(coefs)) 1])));
coefs_zeta = cellfun(zeta, coefs_by_day);
subplot(3,1,2);
% Set axes http://uk.mathworks.com/matlabcentral/newsreader/view_thread/148694
[ax,h1,h2] = plotyy(x,index,x,coefs_zeta);
set(ax(2),'YLim',[0 0.5]);
%%

subplot(3,1,3);
hold off
plot(x, median_ccfs,'r');
hold on
plot(x,tenthP_ccfs,'g');
hold on
plot(x,ninethyP_ccfs,'b');

end