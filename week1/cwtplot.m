%%
figure
subplot(3,1,1)
coefs = cwt(y,1:32,'meyr','plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');


%%

coefs_by_day = num2cell(coefs,1);
morethan = @(tresh) @(x) min(sum(x>tresh)/length(x),0.5);
eta = morethan(1);
coefs_eta = cellfun(eta, coefs_by_day);
subplot(3,1,2);
[ax,h1,h2] = plotyy(x,index,x,coefs_eta);
set(ax(2),'YLim',[0 0.5]);
%%

pfunc = @(q) @(x) prctile(x,q);
ninethyP = pfunc(90);
tenthP = pfunc(10);

tenthP_ccfs = cellfun(tenthP, coefs_by_day);
median_ccfs = cellfun(@median, coefs_by_day);
ninethyP_ccfs = cellfun(ninethyP, coefs_by_day);

subplot(3,1,3);
hold off
plot(x, median_ccfs,'r');
hold on
plot(x,tenthP_ccfs,'g');
hold on
plot(x,ninethyP_ccfs,'b');