%%
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, zetas_mean10);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('Mean 10 days')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(testTS.dates, sign(fts2mat(spreadTS.spread)), testTS.dates, zetas_mean10);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
pair = [zetas_mean10, fts2mat(spreadTS.spread)];
pair(any(isnan(pair),2),:) = [];
crosscorr(pair(:,1),pair(:,2),100);

%%
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('Delay 5 days')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(testTS.dates, sign(fts2mat(spreadTS.spread)), testTS.dates, zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
pair = [zetas_delay5, fts2mat(spreadTS.spread)];
pair(any(isnan(pair),2),:) = [];
crosscorr(pair(:,1),pair(:,2),100);