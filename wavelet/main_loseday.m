indicator_raw = fts2mat(spreadTS.spread) > 5;
indicator_smooth = indicator_raw;

for ii=1:length(spreadTS)
    if sum(indicator_raw(max(ii-5,1):ii-1)) > 0
        indicator_smooth(ii) = 0 ;
    end
end

% losedays = spreadTS.dates(indicator_smooth);
% chk = fts_window(RANGE(indicator_smooth));
% losedays == chk.dates

LOSSRANGE = RANGE(indicator_smooth) - 20; % backwards on indice = trading days

for date = LOSSRANGE
   coefs = calcwt(fts_window(date-windowSize:date));
   figure;
   colormap(pink(240))

   axeAct = subplot(1,1,1);
   image(wcodemat(abs(coefs),240,'mat',1))
   set(axeAct,'YDir','normal')
end

subplot(3,1,1);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, indicator_raw);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[-2,2])
legend([h2;h1],'zeta','spread')


subplot(3,1,2);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, indicator_smooth);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[-2,2])
legend([h2;h1],'zeta','spread')
pair = [zetas_delay5, fts2mat(spreadTS.spread)];
pair(any(isnan(pair),2),:) = [];
crosscorr(pair(:,1),pair(:,2),100);

figure;
subplot(2,1,1);
t = 1;
[hxs,h1,h2]=plotyy(testTS.dates(t:t+500), fts2mat(spreadTS.spread(t:t+500)), testTS.dates(t:t+500), zetas_delay5(t:t+500));
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,1.5])
%set(hxs(),'YLim',[-2,2])
legend([h2;h1],'zeta','spread')

subplot(2,1,2);
t = 1;
zetas_lag = [NaN(8,1);zetas_mean10];
zetas_lag(end-9:end) = [];
[hxs,h1,h2]=plotyy(testTS.dates(t:t+500), fts2mat(spreadTS.spread(t:t+500)), testTS.dates(t:t+500), zetas_lag(t:t+500));
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,1.5])
%set(hxs(),'YLim',[-2,2])
legend([h2;h1],'zeta','spread')


